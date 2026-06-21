-- ============================================================
-- feature_engineering.sql
--
-- Standalone reference copy of the core SQL used in this project.
-- (The same logic is also embedded in src/run_knockout_screen.py,
-- where it's executed via pandas.read_sql_query -- this file exists
-- so the SQL itself is reviewable without reading through Python.)
--
-- Run interactively against data/metabolic_model.db, e.g.:
--   sqlite3 data/metabolic_model.db < sql/feature_engineering.sql
-- ============================================================


-- ------------------------------------------------------------
-- 1. Reaction connectivity: genes, metabolites, reactants, products
-- ------------------------------------------------------------
SELECT
    r.abbreviation,
    r.subsystem,
    COALESCE(gene_counts.n_genes, 0)        AS n_genes,
    COALESCE(met_counts.n_metabolites, 0)   AS n_metabolites,
    COALESCE(reactant_counts.n_reactants, 0) AS n_reactants,
    COALESCE(product_counts.n_products, 0)   AS n_products
FROM reactions r
LEFT JOIN (
    SELECT reaction_abbr, COUNT(DISTINCT gene_id) AS n_genes
    FROM reaction_genes
    GROUP BY reaction_abbr
) gene_counts ON gene_counts.reaction_abbr = r.abbreviation
LEFT JOIN (
    SELECT reaction_abbr, COUNT(DISTINCT metabolite_abbr) AS n_metabolites
    FROM reaction_metabolites
    GROUP BY reaction_abbr
) met_counts ON met_counts.reaction_abbr = r.abbreviation
LEFT JOIN (
    SELECT reaction_abbr, COUNT(*) AS n_reactants
    FROM reaction_metabolites WHERE role = 'reactant'
    GROUP BY reaction_abbr
) reactant_counts ON reactant_counts.reaction_abbr = r.abbreviation
LEFT JOIN (
    SELECT reaction_abbr, COUNT(*) AS n_products
    FROM reaction_metabolites WHERE role = 'product'
    GROUP BY reaction_abbr
) product_counts ON product_counts.reaction_abbr = r.abbreviation
ORDER BY n_genes DESC
LIMIT 20;


-- ------------------------------------------------------------
-- 2. Subsystems ranked by reaction count
-- ------------------------------------------------------------
SELECT subsystem, COUNT(*) AS n_reactions
FROM reactions
GROUP BY subsystem
ORDER BY n_reactions DESC;


-- ------------------------------------------------------------
-- 3. Genes that participate in the most reactions
--    (candidate hub genes / pleiotropy proxy)
-- ------------------------------------------------------------
SELECT gene_id, COUNT(*) AS n_reactions
FROM reaction_genes
GROUP BY gene_id
ORDER BY n_reactions DESC
LIMIT 20;


-- ------------------------------------------------------------
-- 4. Metabolites that participate in the most reactions
--    (hub metabolites -- e.g. ATP, water, NAD+ typically top this list)
-- ------------------------------------------------------------
SELECT
    m.abbreviation,
    m.description,
    COUNT(*) AS n_reactions
FROM reaction_metabolites rm
JOIN metabolites m ON m.abbreviation = rm.metabolite_abbr
GROUP BY m.abbreviation, m.description
ORDER BY n_reactions DESC
LIMIT 20;


-- ------------------------------------------------------------
-- 5. Reactions per compartment-pair (transport reaction summary)
-- ------------------------------------------------------------
SELECT
    is_transport,
    is_exchange,
    is_reversible,
    COUNT(*) AS n_reactions
FROM reactions
GROUP BY is_transport, is_exchange, is_reversible
ORDER BY n_reactions DESC;


-- ------------------------------------------------------------
-- 6. Window function example: rank reactions within each subsystem
--    by number of genes (useful for spotting the "most regulated"
--    reaction in every pathway)
-- ------------------------------------------------------------
SELECT
    abbreviation,
    subsystem,
    n_genes,
    RANK() OVER (
        PARTITION BY subsystem
        ORDER BY n_genes DESC
    ) AS rank_within_subsystem
FROM (
    SELECT
        r.abbreviation,
        r.subsystem,
        COUNT(DISTINCT rg.gene_id) AS n_genes
    FROM reactions r
    LEFT JOIN reaction_genes rg ON rg.reaction_abbr = r.abbreviation
    GROUP BY r.abbreviation, r.subsystem
)
WHERE n_genes > 0
ORDER BY subsystem, rank_within_subsystem
LIMIT 30;
