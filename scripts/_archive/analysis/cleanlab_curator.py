"""
===============================================================================
Module: cleanlab_curator.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Confident Learning algorithms (Cleanlab) and cross-validated out-of-fold
    probability estimators to identify taxonomic misidentification noise.
===============================================================================
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold

logger = logging.getLogger("CleanlabCurator")


def compute_out_of_fold_probabilities(
    features: np.ndarray,
    labels: np.ndarray,
    n_splits: int = 5,
    random_state: int = 42
) -> Tuple[np.ndarray, float, float]:
    """
    Fits stratified cross-validated Logistic Regression models to compute out-of-fold probabilities.
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    num_classes = len(np.unique(labels))
    pred_probs = np.zeros((len(labels), num_classes), dtype=np.float64)

    y_true_all, y_pred_all = [], []

    for train_idx, val_idx in skf.split(features, labels):
        X_train, y_train = features[train_idx], labels[train_idx]
        X_val, y_val = features[val_idx], labels[val_idx]

        clf = LogisticRegression(max_iter=1000, C=1.0, random_state=random_state)
        clf.fit(X_train, y_train)

        probs = clf.predict_proba(X_val)
        pred_probs[val_idx] = probs

        y_true_all.extend(y_val)
        y_pred_all.extend(np.argmax(probs, axis=1))

    acc = accuracy_score(y_true_all, y_pred_all)
    f1 = f1_score(y_true_all, y_pred_all, average="weighted")
    logger.info(f"OOF Classifier Cross-Validation: Accuracy = {acc:.4f}, Weighted F1 = {f1:.4f}")

    return pred_probs, acc, f1


def run_confident_learning_audit(
    pred_probs: np.ndarray,
    labels: np.ndarray,
    records_df: pd.DataFrame,
    class_names: List[str],
    error_threshold: float = 0.85
) -> pd.DataFrame:
    """
    Applies Confident Learning algorithms to identify potential herbarium label noise.
    """
    try:
        import cleanlab
        from cleanlab.filter import find_label_issues
        from cleanlab.rank import get_label_quality_scores

        quality_scores = get_label_quality_scores(labels=labels, pred_probs=pred_probs)
        issues = find_label_issues(labels=labels, pred_probs=pred_probs, return_indices_ranked_by="self_confidence")
    except ImportError:
        logger.warning("Cleanlab library not available; computing heuristic label error margins.")
        predicted_classes = np.argmax(pred_probs, axis=1)
        quality_scores = np.array([pred_probs[i, labels[i]] for i in range(len(labels))])
        issues = np.where((predicted_classes != labels) & (1.0 - quality_scores > error_threshold))[0]

    audit_df = records_df.copy()
    num_records = len(audit_df)

    given_classes = [class_names[l] if l < len(class_names) else "Unknown" for l in labels]
    pred_indices = np.argmax(pred_probs, axis=1)
    pred_classes = [class_names[i] if i < len(class_names) else "Unknown" for i in pred_indices]

    conf_given = np.array([pred_probs[i, labels[i]] if i < len(pred_probs) and labels[i] < pred_probs.shape[1] else 0.0 for i in range(num_records)])
    conf_pred = np.max(pred_probs, axis=1)
    c_error = 1.0 - conf_given

    audit_df["species_raw"] = audit_df.get("species_raw", audit_df.get("species", given_classes))
    audit_df["species_standardized"] = given_classes
    audit_df["given_label"] = given_classes
    audit_df["predicted_label"] = pred_classes
    audit_df["confidence_given_class"] = np.round(conf_given, 4)
    audit_df["confidence_predicted_class"] = np.round(conf_pred, 4)
    audit_df["label_quality_score"] = np.round(quality_scores, 4)
    audit_df["c_error"] = np.round(c_error, 4)

    is_issue = np.zeros(num_records, dtype=bool)
    is_issue[issues] = True
    audit_df["is_cleanlab_issue"] = is_issue
    audit_df["is_label_corrupted"] = audit_df["c_error"] > error_threshold

    triage_actions = []
    reasons = []
    for idx, row in audit_df.iterrows():
        if row["is_label_corrupted"]:
            triage_actions.append("Prune & Queue for Annotation Triage")
            reasons.append(f"Deep vision DINOv2 predicts {row['predicted_label']} (conf: {row['confidence_predicted_class']:.2f}) vs recorded {row['given_label']}")
        elif row["is_cleanlab_issue"]:
            triage_actions.append("Flag for Morphometric Review")
            reasons.append(f"Confident learning flags potential discordance (quality: {row['label_quality_score']:.2f})")
        else:
            triage_actions.append("Retain")
            reasons.append("High label consistency across visual self-supervised embeddings")

    audit_df["triage_action"] = triage_actions
    audit_df["discordance_reason"] = reasons

    flagged_count = int(audit_df["is_label_corrupted"].sum())
    logger.info(f"Confident Learning Audit completed: Flagged {flagged_count} / {len(audit_df)} corrupted labels.")
    return audit_df
