"""
export_cobra_json.py
----------------------
Exports the SQLite-stored model as a standard COBRA model JSON file
(the same format produced by cobrapy's `cobra.io.save_json_model`).

Why this matters for Escher:
    This Drosophila model has no pre-built Escher map (unlike e.g.
    E. coli core metabolism or Human-GEM, which ship official maps).
    To visualize it in Escher you have two options:

      1. Open Escher Builder -> Map > "New map" -> then
         Model > Load COBRA model JSON (this file) -> Escher will let
         you build a map by dragging reactions in from the loaded
         model, auto-laying out the ones you pick.

      2. Use Escher's "Build new map from model" workflow (in the
         Escher Python package, escher.Builder(model_json=...)), which
         can auto-generate a starter layout from the whole model or a
         chosen subsystem -- recommended given this model has 8,000+
         reactions, so picking ONE subsystem at a time (e.g. just
         "Glycolysis / Gluconeogenesis") will give a far more readable
         map than the entire network at once.

Output: output/drosophila_cobra_model.json
"""

import json
import sqlite3
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
OUTPUT_DIR = HERE.parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

DB_PATH = DATA_DIR / "metabolic_model.db"
OUT_PATH = OUTPUT_DIR / "drosophila_cobra_model.json"


def main():
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    # --- metabolites ---
    metabolites = []
    for abbr, desc, formula, compartment in c.execute(
        "SELECT abbreviation, description, formula, compartment FROM metabolites"
    ):
        # Escher/COBRA convention: id ends in _<compartment letter>, e.g. "glc__D_c"
        compartment_letter = abbr.split("[")[-1].rstrip("]") if "[" in abbr else "c"
        cobra_id = abbr.replace("[", "_").replace("]", "")
        metabolites.append(
            {
                "id": cobra_id,
                "name": desc or abbr,
                "compartment": compartment_letter,
                "formula": formula or "",
            }
        )

    # --- reactions ---
    reactions = []
    rxn_rows = list(
        c.execute(
            "SELECT abbreviation, description, lower_bound, upper_bound, subsystem FROM reactions"
        )
    )
    for abbr, desc, lb, ub, subsystem in rxn_rows:
        # metabolite stoichiometry for this reaction
        stoich = {}
        for met_abbr, coeff in c.execute(
            "SELECT metabolite_abbr, stoichiometry FROM reaction_metabolites WHERE reaction_abbr = ?",
            (abbr,),
        ):
            cobra_met_id = met_abbr.replace("[", "_").replace("]", "")
            stoich[cobra_met_id] = coeff

        # gene reaction rule, rebuilt as "(gene1 or gene2)" -- OR-only
        # approximation; the original AND/OR boolean structure isn't
        # retained at the SQL layer (reaction_genes is a flat many-to-many
        # table), so this is a simplification good enough for map display
        # but not for strict GPR logic evaluation.
        genes = [
            row[0]
            for row in c.execute(
                "SELECT gene_id FROM reaction_genes WHERE reaction_abbr = ?",
                (abbr,),
            )
        ]
        gpr = " or ".join(genes) if genes else ""

        reactions.append(
            {
                "id": abbr,
                "name": desc or abbr,
                "metabolites": stoich,
                "lower_bound": lb,
                "upper_bound": ub,
                "gene_reaction_rule": gpr,
                "subsystem": subsystem or "",
            }
        )

    # --- genes ---
    genes_list = [
        {"id": row[0], "name": row[0]}
        for row in c.execute("SELECT fly_base_id FROM genes")
    ]

    model = {
        "id": "Drosophila_metabolic_model",
        "compartments": {
            "c": "Cytosol",
            "m": "Mitochondria",
            "r": "Endoplasmic reticulum",
            "e": "Extracellular",
            "p": "Peroxisome",
            "n": "Nucleus",
            "l": "Lysosome",
            "g": "Golgi apparatus",
            "i": "Mitochondrial intermembrane space",
            "s": "Extracellular",  # exchange-boundary pseudo-compartment seen in this file
        },
        "metabolites": metabolites,
        "reactions": reactions,
        "genes": genes_list,
        "version": "1",
    }

    with open(OUT_PATH, "w") as f:
        json.dump(model, f)

    print(f"Wrote COBRA-format model JSON: {OUT_PATH}")
    print(f"  {len(reactions)} reactions, {len(metabolites)} metabolites, {len(genes_list)} genes")
    print()
    print("NOTE: reaction IDs in the data files (escher_*.json) use the raw")
    print("abbreviation (e.g. 'HMR_3905'); metabolite IDs in this model JSON")
    print("use underscores instead of brackets (e.g. 'm02552_c') to match")
    print("standard COBRA/Escher metabolite ID conventions.")

    conn.close()


if __name__ == "__main__":
    main()
