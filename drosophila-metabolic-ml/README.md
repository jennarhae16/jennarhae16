# Predicting Metabolic Reaction Essentiality in a *Drosophila* Genome-Scale Model

This project takes a real genome-scale metabolic reconstruction for
*Drosophila melanogaster* (8,230 reactions, 6,990 metabolites, 2,388
genes) and asks: **can cheap, SQL-derived structural features predict
which reactions are essential for growth, without ever running flux
balance analysis (FBA) on them?**

The pipeline goes: raw spreadsheet model &rarr; normalized SQL database
&rarr; SQL feature engineering &rarr; FBA-derived ground-truth labels
&rarr; ML classifier &rarr; Escher-compatible visualization export.

## Why this is interesting

Genome-scale metabolic models encode the full reaction network of an
organism's metabolism, including which genes catalyze which reactions
(gene-protein-reaction rules). A standard question in systems biology
is *reaction/gene essentiality*: if you knock out a reaction, does the
organism's simulated growth rate collapse? Answering this directly
requires solving a linear program (FBA) for every single reaction.

This project tests whether you can shortcut that: using only
network-topology features you can compute with plain SQL (gene count,
metabolite connectivity, subsystem, reversibility, transport/exchange
flags) — features that say nothing about flux *values*, only about
*structure* — how well can a classifier flag the reactions that FBA
says are essential?

## Pipeline

```
data/Drosophila.xlsx                 (raw model: Reaction List + Metabolite List sheets)
        |
        v
code/build_database.py                parses equations & GPR rules -> normalized SQLite DB
        |
        v
data/metabolic_model.db              5 tables: reactions, metabolites, genes,
        |                             reaction_genes, reaction_metabolites
        v
code/fba.py                           FBA engine (scipy.optimize.linprog, sparse S matrix)
        |
        v
code/run_knockout_screen.py           SQL feature queries + full single-reaction
        |                             knockout screen -> ground-truth essential/non-essential
        v
data/reaction_features_labeled.csv
        |
        v
code/train_model.py                   Random Forest classifier, SQL features -> essentiality
        |
        v
output/model_performance.json, feature_importance.csv, trained_model.pkl
        |
        v
code/export_escher_data.py            -> Escher-ready reaction data JSON files
code/export_cobra_json.py             -> COBRA model JSON (for building a custom Escher map)
```

## Database schema

```
metabolites          (abbreviation, description, formula, compartment)
reactions            (abbreviation, description, equation, lower_bound,
                       upper_bound, is_objective, subsystem)
genes                (fly_base_id)
reaction_genes       (reaction_abbr, gene_id)            -- GPR rule
reaction_metabolites (reaction_abbr, metabolite_abbr,
                       stoichiometry, role)              -- one row per term
```


## Essentially Discovery by Knockout

For each reaction, its flux is forced to zero and FBA is re-solved.
A reaction is labeled **essential** if optimal growth drops by more
than 1% relative to wild type.

**Speed note:** any reaction carrying zero flux at the wild-type
optimum cannot affect the objective when additionally forced to zero,
so those ~6,400 reactions are labeled non-essential directly and only
the ~1,800 active reactions are actually re-solved through the LP.
This is exact, not an approximation — it just skips LP calls that are
guaranteed to return the wild-type value.

**Result:** 258 of 8,229 non-biomass reactions (3.1%) are essential —
consistent with the high redundancy typical of genome-scale metabolic
networks, where most reactions have backup routes.

## ML model

A Random Forest classifier (`class_weight="balanced"` to handle the
~3% positive rate) trained on:

- gene count, metabolite count, reactant/product count
- reversibility, exchange flag, transport flag
- subsystem (one-hot, top 20 + "Other")
- flux bound range

**Test set performance** (25% holdout, stratified):

| metric                | value |
|------------------------|-------|
| ROC AUC                | 0.857 |
| Recall (essential)     | 0.58  |
| Precision (essential)  | 0.17  |
| Accuracy               | 0.89  |

The low precision is expected and reported honestly: essential
reactions are rare (~3%), so even a useful model produces many false
positives in absolute terms. The meaningful number here is recall —
the model flags over half of the truly essential reactions using only
structural features, with zero FBA solves needed at inference time.

Top predictive features: subsystem membership (Glycerolipid
metabolism in particular), gene count, metabolite connectivity, and
flux bound range. See `output/feature_importance.csv` for the full
ranking.

## Escher visualization export

`src/export_escher_data.py` produces three Escher-compatible reaction
data files (the simple `{reaction_id: value}` JSON format Escher's
**Data &rarr; Load reaction data** menu accepts directly):

| file                                    | values                                          |
|-------------------------------------------|--------------------------------------------------|
| `escher_true_essentiality.json`           | 1.0 / 0.0 — FBA ground truth                      |
| `escher_predicted_probability.json`       | 0–1 — model's predicted essentiality probability  |
| `escher_growth_ratio.json`                | 0–1 — fraction of WT growth retained post-knockout |

**This model has no pre-built Escher map** (unlike *E. coli* core
metabolism or Human-GEM, which ship official ones), so to visualize
these values a map needed to be built.
![Glycolysis Pathway](https://github.com/jennarhae16/jennarhae16/blob/main/drosophila-metabolic-ml/output/Glycolysis_Dros.png)

Note: metabolite IDs in `drosophila_cobra_model.json` use underscores
instead of brackets (`m02552_c` instead of `m02552[c]`) to match
standard COBRA/Escher ID conventions; reaction IDs are unchanged.

## Limitations / honest caveats

- The gene-reaction rule (GPR) boolean logic (AND/OR) is simplified
  when exploded into the `reaction_genes` table — the relational
  schema captures *which* genes are associated with a reaction but
  not the AND/OR structure, so the COBRA JSON export reconstructs GPR
  strings as OR-only. This does not affect the FBA knockout screen
  (which operates on reactions directly) but would need fixing for
  rigorous single-gene-deletion analysis.
- This is a single organism, single growth condition. The "essential"
  label is specific to the biomass objective and default exchange
  bounds encoded in the original model — change the medium and some
  labels would change too.
- Low precision (17%) on the positive class means this model is
  useful as a triage/prioritization tool (which reactions deserve a
  real FBA check first), not a replacement for running FBA itself.
