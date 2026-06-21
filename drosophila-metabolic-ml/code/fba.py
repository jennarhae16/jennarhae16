"""
fba.py
------
A minimal Flux Balance Analysis (FBA) engine built directly on
scipy.optimize.linprog, reading the reaction network from the SQLite
database produced by build_database.py.

This avoids a dependency on cobrapy (not installable in this sandbox)
while implementing the same underlying LP:

    maximize     c^T v
    subject to   S v = 0          (steady state: production = consumption)
                 lb <= v <= ub    (flux bounds)

where v is the vector of reaction fluxes, S is the stoichiometric
matrix (metabolites x reactions), and c is the objective vector
(1 for the biomass reaction, 0 elsewhere).

scipy.optimize.linprog only minimizes, so we minimize -c^T v.
"""

import sqlite3
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
DB_PATH = DATA_DIR / "metabolic_model.db"


class MetabolicModel:
    """Loads the reaction network from SQLite and builds the LP matrices."""

    def __init__(self, db_path: Path = DB_PATH):
        self.conn = sqlite3.connect(db_path)
        self._load()

    def _load(self):
        c = self.conn.cursor()

        self.reactions = [
            row[0]
            for row in c.execute(
                "SELECT abbreviation FROM reactions ORDER BY abbreviation"
            )
        ]
        self.rxn_index = {r: i for i, r in enumerate(self.reactions)}

        self.metabolites = [
            row[0]
            for row in c.execute(
                "SELECT abbreviation FROM metabolites ORDER BY abbreviation"
            )
        ]
        self.met_index = {m: i for i, m in enumerate(self.metabolites)}

        n_rxn = len(self.reactions)
        n_met = len(self.metabolites)

        self.lower_bound = np.zeros(n_rxn)
        self.upper_bound = np.zeros(n_rxn)
        self.objective = np.zeros(n_rxn)
        self.subsystem = [None] * n_rxn

        for abbr, lb, ub, is_obj, subsystem in c.execute(
            "SELECT abbreviation, lower_bound, upper_bound, is_objective, subsystem FROM reactions"
        ):
            i = self.rxn_index[abbr]
            self.lower_bound[i] = lb
            self.upper_bound[i] = ub
            self.objective[i] = 1.0 if is_obj else 0.0
            self.subsystem[i] = subsystem

        # Sparse stoichiometric matrix S (metabolites x reactions)
        S = lil_matrix((n_met, n_rxn))
        for rxn_abbr, met_abbr, coeff in c.execute(
            "SELECT reaction_abbr, metabolite_abbr, stoichiometry FROM reaction_metabolites"
        ):
            if met_abbr not in self.met_index or rxn_abbr not in self.rxn_index:
                continue
            r = self.met_index[met_abbr]
            col = self.rxn_index[rxn_abbr]
            S[r, col] += coeff
        self.S = S.tocsr()

        # Map reaction -> list of genes (for single-gene-deletion screens)
        self.reaction_genes = {r: [] for r in self.reactions}
        for rxn_abbr, gene_id in c.execute(
            "SELECT reaction_abbr, gene_id FROM reaction_genes"
        ):
            if rxn_abbr in self.reaction_genes:
                self.reaction_genes[rxn_abbr].append(gene_id)

    def solve(self, lower_bound=None, upper_bound=None):
        """Solve the FBA LP. Returns (objective_value, flux_vector) or (0.0, None) if infeasible.

        lower_bound / upper_bound: optional override arrays (same length as
        self.reactions) used for knockout simulations.
        """
        lb = self.lower_bound if lower_bound is None else lower_bound
        ub = self.upper_bound if upper_bound is None else upper_bound

        n_rxn = len(self.reactions)
        n_met = len(self.metabolites)

        # linprog minimizes; we want to maximize objective^T v -> minimize -objective^T v
        cost = -self.objective

        # Equality constraint: S v = 0
        A_eq = self.S
        b_eq = np.zeros(n_met)

        bounds = list(zip(lb, ub))

        result = linprog(
            cost,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
        )

        if not result.success:
            return 0.0, None

        objective_value = -result.fun
        return objective_value, result.x

    def knockout_reaction(self, reaction_abbr):
        """Simulate a single-reaction knockout: force its flux to 0."""
        i = self.rxn_index[reaction_abbr]
        lb = self.lower_bound.copy()
        ub = self.upper_bound.copy()
        lb[i] = 0.0
        ub[i] = 0.0
        return self.solve(lower_bound=lb, upper_bound=ub)

    def knockout_gene(self, gene_id):
        """Simulate a single-gene knockout: zero out every reaction whose GPR
        depends solely on this gene (simple OR/AND handling: a reaction is
        knocked out only if ALL of its genes are removed, i.e. we approximate
        AND-linked complexes needing every subunit, and OR-linked isozymes
        as redundant unless this is the only gene).

        For simplicity (and because most GPRs in this model are single-gene
        or pure-OR), we treat: reaction lost if gene_id is in its gene list
        AND it is the reaction's only annotated gene OR all of the reaction's
        genes are this same gene.
        """
        affected = [
            r
            for r, genes in self.reaction_genes.items()
            if genes == [gene_id]  # reaction depends only on this single gene
        ]
        lb = self.lower_bound.copy()
        ub = self.upper_bound.copy()
        for r in affected:
            i = self.rxn_index[r]
            lb[i] = 0.0
            ub[i] = 0.0
        obj, flux = self.solve(lower_bound=lb, upper_bound=ub)
        return obj, flux, affected


if __name__ == "__main__":
    model = MetabolicModel()
    print(f"Loaded model: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites")

    wt_growth, wt_flux = model.solve()
    print(f"Wild-type optimal growth (biomass flux): {wt_growth:.6f}")

    if wt_flux is not None:
        active = np.sum(np.abs(wt_flux) > 1e-6)
        print(f"Active reactions at optimum: {active} / {len(model.reactions)}")
