"""
export_escher_data.py
-----------------------
Exports model results as Escher "reaction data" JSON files -- the
simple {reaction_id: value} dictionary format that Escher's
Data > Load reaction data menu accepts directly.

This produces THREE separate JSON files so you can load whichever
view you want onto an Escher map of this Drosophila model (or any
BiGG-style map that shares reaction IDs):

  1. escher_true_essentiality.json
        1.0 if FBA-confirmed essential, 0.0 otherwise.
        Use this to highlight the ground-truth essential reactions.

  2. escher_predicted_probability.json
        The ML model's predicted probability of essentiality (0-1),
        for every reaction -- including the ones never explicitly
        knocked out. Use this to see where the model is confident.

  3. escher_growth_ratio.json
        Fraction of wild-type growth retained after knockout (1.0 =
        no effect, 0.0 = lethal). This is the closest analog to a
        flux map and gives the smoothest gradient for coloring.

Escher expects reaction IDs as keys. This model's reaction
abbreviations (e.g. "HMR_3905") match the IDs used in the
Human-GEM-style maps that this Drosophila model was adapted from, so
they should line up directly if you load an HMR/Human1-based Escher
map. If you build a custom map from this model's own reactions
instead, the IDs will already match by construction.
"""

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
OUTPUT_DIR = HERE.parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

PREDICTIONS_PATH = OUTPUT_DIR / "reactions_with_predictions.csv"


def main():
    df = pd.read_csv(PREDICTIONS_PATH)

    # 1. Ground-truth essentiality (from the FBA knockout screen)
    true_essential = dict(zip(df["abbreviation"], df["essential"].astype(float)))
    with open(OUTPUT_DIR / "escher_true_essentiality.json", "w") as f:
        json.dump(true_essential, f)

    # 2. ML-predicted probability of essentiality
    predicted_proba = dict(
        zip(df["abbreviation"], df["predicted_essential_proba"].round(4))
    )
    with open(OUTPUT_DIR / "escher_predicted_probability.json", "w") as f:
        json.dump(predicted_proba, f)

    # 3. Growth ratio after knockout (smooth gradient, good for flux-style coloring)
    growth_ratio = dict(zip(df["abbreviation"], df["growth_ratio"].round(4)))
    with open(OUTPUT_DIR / "escher_growth_ratio.json", "w") as f:
        json.dump(growth_ratio, f)

    print("Escher-compatible JSON files written to output/:")
    print(f"  escher_true_essentiality.json     ({len(true_essential)} reactions)")
    print(f"  escher_predicted_probability.json ({len(predicted_proba)} reactions)")
    print(f"  escher_growth_ratio.json          ({len(growth_ratio)} reactions)")
    print()
    print("In Escher: Data menu -> Load reaction data -> pick one of these files.")


if __name__ == "__main__":
    main()
