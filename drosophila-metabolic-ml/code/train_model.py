"""
train_model.py
---------------
Trains a classifier to predict reaction essentiality (the FBA-derived
ground truth) using only cheap, SQL-computed structural features --
i.e. without ever running FBA on the reaction being predicted.

This tests whether network topology and annotation features (gene
count, connectivity, subsystem, reversibility, exchange/transport
flags) carry enough signal to flag likely-essential reactions without
the cost of solving an LP for every one.

Input : data/reaction_features_labeled.csv
Output: output/model_performance.json
        output/feature_importance.csv
        output/trained_model.pkl
"""

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
OUTPUT_DIR = HERE.parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

IN_PATH = DATA_DIR / "reaction_features_labeled.csv"

NUMERIC_FEATURES = [
    "is_reversible",
    "is_exchange",
    "is_transport",
    "bound_range",
    "n_genes",
    "n_metabolites",
    "n_reactants",
    "n_products",
]
CATEGORICAL_FEATURES = ["subsystem"]
TARGET = "essential"


def load_data():
    df = pd.read_csv(IN_PATH)
    return df


def build_feature_matrix(df: pd.DataFrame, encoder: OneHotEncoder = None, fit=True):
    numeric = df[NUMERIC_FEATURES].fillna(0).values

    # Collapse rare subsystems so one-hot encoding doesn't explode into
    # 126 near-empty columns (most reactions live in fewer than 20 systems
    # that matter for the model).
    top_subsystems = df["subsystem"].value_counts().head(20).index
    subsystem_clipped = df["subsystem"].where(
        df["subsystem"].isin(top_subsystems), other="Other"
    ).values.reshape(-1, 1)

    if fit:
        encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
        cat_encoded = encoder.fit_transform(subsystem_clipped)
    else:
        cat_encoded = encoder.transform(subsystem_clipped)

    X = np.hstack([numeric, cat_encoded])
    feature_names = NUMERIC_FEATURES + list(encoder.get_feature_names_out(["subsystem"]))
    return X, feature_names, encoder


def main():
    df = load_data()
    y = df[TARGET].values

    X, feature_names, encoder = build_feature_matrix(df, fit=True)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=42, stratify=y
    )

    print(f"Train size: {len(y_train)}  Test size: {len(y_test)}")
    print(f"Essential rate -- train: {y_train.mean():.3f}  test: {y_test.mean():.3f}")

    model = RandomForestClassifier(
        n_estimators=300,
        max_depth=12,
        class_weight="balanced",  # essential reactions are rare (~3%)
        random_state=42,
        n_jobs=-1,
    )
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    report = classification_report(y_test, y_pred, output_dict=True)
    auc = roc_auc_score(y_test, y_proba)
    cm = confusion_matrix(y_test, y_pred).tolist()

    print()
    print(classification_report(y_test, y_pred))
    print(f"ROC AUC: {auc:.3f}")
    print(f"Confusion matrix: {cm}")

    # Feature importance
    importance = pd.DataFrame(
        {"feature": feature_names, "importance": model.feature_importances_}
    ).sort_values("importance", ascending=False)

    print()
    print("Top 10 most important features:")
    print(importance.head(10).to_string(index=False))

    # Save outputs
    with open(OUTPUT_DIR / "model_performance.json", "w") as f:
        json.dump(
            {
                "roc_auc": auc,
                "confusion_matrix": cm,
                "classification_report": report,
                "train_size": len(y_train),
                "test_size": len(y_test),
            },
            f,
            indent=2,
        )

    importance.to_csv(OUTPUT_DIR / "feature_importance.csv", index=False)

    with open(OUTPUT_DIR / "trained_model.pkl", "wb") as f:
        pickle.dump({"model": model, "encoder": encoder, "feature_names": feature_names}, f)

    # Also save full predictions (all reactions, not just test set) for the
    # Escher export step -- predicted probability of essentiality per
    # reaction, plus the true FBA label, side by side.
    X_all, _, _ = build_feature_matrix(df, encoder=encoder, fit=False)
    df["predicted_essential_proba"] = model.predict_proba(X_all)[:, 1]
    df.to_csv(OUTPUT_DIR / "reactions_with_predictions.csv", index=False)

    print(f"\nSaved model, metrics, and predictions to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
