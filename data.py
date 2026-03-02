from __future__ import annotations

import os
import re
from dataclasses import dataclass
from io import StringIO
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import requests
import streamlit as st

from config import DATASETS


ANTIGEN_CANDIDATE_COLUMNS = [
    # Generic antigen / gene / symbol names
    "antigen",
    "Antigen",
    "gene",
    "Gene",
    "Gene_Names",
    "gene_name",
    "GeneName",
    "symbol",
    "Symbol",
    "target",
    "Target",
    "ID",
    "id",
    # Common feature identifiers in these datasets
    "Feature_ID",
    "feature_id",
    "Accession",
    "accession",
    # Descriptive protein/antigen names
    "Protein_Names",
    "protein_name",
    "Protein",
    "protein",
    "Name",
    "name",
    # Sometimes species name is the antigen-like unit
    "Species",
    "species",
]

LOGFC_CANDIDATE_COLUMNS = ["logFC", "logfc", "LogFC"]
FDR_CANDIDATE_COLUMNS = ["FDR", "fdr", "adj.P.Val", "adj.P.Val.", "padj", "qvalue", "q_value"]
PVALUE_CANDIDATE_COLUMNS = ["PValue", "pvalue", "P.Value", "P.Value.", "p_val", "pvalue"]
LOGCPM_CANDIDATE_COLUMNS = ["logCPM", "logcpm", "LogCPM"]
AVE_EXPR_CANDIDATE_COLUMNS = ["AveExpr", "ave_expr", "Amean", "avg_log_expr"]


@dataclass
class ResolvedSource:
    source_type: str  # "path" or "url"
    value: str
    description: str


def _normalize_key(s: str) -> str:
    """Uppercase and replace non-alphanumeric with underscores for ENV var compatibility."""
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").upper()


def _env_var_candidates(dataset_key: str, contrast: Optional[str] = None) -> Dict[str, str]:
    """
    Build possible environment variable names for a dataset (and optional contrast).

    Priority order (checked top to bottom):
    1. DATASET_CONTRAST_CSV_URL
    2. DATASET_CONTRAST_CSV_PATH
    3. DATASET_CSV_URL
    4. DATASET_CSV_PATH
    """
    base = _normalize_key(dataset_key)
    envs: Dict[str, str] = {}
    if contrast:
        c_part = _normalize_key(contrast)
        envs[f"{base}_{c_part}_CSV_URL"] = "url"
        envs[f"{base}_{c_part}_CSV_PATH"] = "path"
    envs[f"{base}_CSV_URL"] = "url"
    envs[f"{base}_CSV_PATH"] = "path"
    return envs


def resolve_contrast_source(dataset_key: str, contrast: str) -> Tuple[Optional[ResolvedSource], str]:
    """
    Resolve a CSV source for a given dataset+contrast using:
    1. Environment variable overrides (with contrast-specific names preferred)
    2. config.DATASETS entries for that dataset+contrast
    """
    attempted: list[str] = []

    # 1) Check env vars
    env_candidates = _env_var_candidates(dataset_key, contrast)
    for env_name, kind in env_candidates.items():
        val = os.getenv(env_name)
        if val:
            src = ResolvedSource(
                source_type="url" if kind == "url" else "path",
                value=val,
                description=f"ENV[{env_name}]",
            )
            return src, ""
        attempted.append(env_name)

    # 2) Fallback to config
    ds_cfg = DATASETS.get(dataset_key, {})
    c_cfg = ds_cfg.get("contrasts", {}).get(contrast, {})
    csv_path = c_cfg.get("csv")
    csv_url = c_cfg.get("csv_url")
    if csv_url:
        return (
            ResolvedSource(source_type="url", value=csv_url, description=f"config.csv_url ({dataset_key}/{contrast})"),
            "",
        )
    if csv_path:
        return (
            ResolvedSource(source_type="path", value=csv_path, description=f"config.csv ({dataset_key}/{contrast})"),
            "",
        )

    msg = "Unable to resolve CSV source. Tried environment variables: " + ", ".join(attempted)
    return None, msg


@st.cache_data(show_spinner=False)
def _read_csv_from_source(source: ResolvedSource) -> pd.DataFrame:
    if source.source_type == "path":
        return pd.read_csv(source.value)
    if source.source_type == "url":
        resp = requests.get(source.value)
        resp.raise_for_status()
        return pd.read_csv(StringIO(resp.text))
    raise ValueError(f"Unsupported source_type: {source.source_type}")


def _pick_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for name in candidates:
        if name in df.columns:
            return name
    # Case-insensitive pass
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def _handle_duplicates(df: pd.DataFrame, strategy: str) -> pd.DataFrame:
    """
    Handle duplicate antigens.

    strategy:
      - "min_fdr": keep row with lowest FDR (default)
      - "max_abs_logfc": keep row with largest |logFC|
      - "mean": aggregate numeric columns as mean
    """
    if "antigen" not in df.columns:
        return df
    if strategy == "mean":
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        non_numeric_cols = [c for c in df.columns if c not in numeric_cols and c != "antigen"]
        agg_spec = {col: "mean" for col in numeric_cols}
        result = df.groupby("antigen", as_index=False).agg(agg_spec)
        # Preserve a representative value for non-numeric columns where possible
        for col in non_numeric_cols:
            first_vals = df.groupby("antigen")[col].first()
            result[col] = result["antigen"].map(first_vals)
        return result

    if strategy == "max_abs_logfc" and "logFC" in df.columns:
        df = df.copy()
        df["_abs_logFC"] = df["logFC"].abs()
        df = df.sort_values(["antigen", "_abs_logFC"], ascending=[True, False])
        result = df.drop_duplicates(subset=["antigen"], keep="first").drop(columns=["_abs_logFC"])
        return result

    # Default: min_fdr
    if "FDR" in df.columns:
        df = df.copy()
        df = df.sort_values(["antigen", "FDR"], ascending=[True, True])
        result = df.drop_duplicates(subset=["antigen"], keep="first")
        return result

    return df.drop_duplicates(subset=["antigen"])


@st.cache_data(show_spinner=False)
def load_and_standardize(
    dataset_key: str,
    contrast: str,
    duplicate_strategy: str = "min_fdr",
) -> Tuple[Optional[pd.DataFrame], Optional[str]]:
    """
    Load and standardize an edgeR results CSV for a given dataset+contrast.

    Returns (DataFrame, error_message). On success, error_message is None.
    """
    source, err = resolve_contrast_source(dataset_key, contrast)
    if source is None:
        return None, err

    try:
        raw_df = _read_csv_from_source(source)
    except Exception as exc:  # noqa: BLE001
        return None, f"Failed to load CSV from {source.description}: {exc}"

    if raw_df.empty:
        return None, "Loaded CSV is empty."

    df = raw_df.copy()

    # Antigen column
    antigen_col = _pick_column(df, ANTIGEN_CANDIDATE_COLUMNS)
    if antigen_col is None:
        return None, "Could not find an antigen column in CSV."
    df = df.rename(columns={antigen_col: "antigen"})

    # Required logFC
    logfc_col = _pick_column(df, LOGFC_CANDIDATE_COLUMNS)
    if logfc_col is None:
        return None, "Could not find a logFC column in CSV."
    if logfc_col != "logFC":
        df = df.rename(columns={logfc_col: "logFC"})

    # Required FDR (or adj.P.Val fallback)
    fdr_col = _pick_column(df, FDR_CANDIDATE_COLUMNS)
    if fdr_col is None:
        return None, "Could not find an FDR/adj.P.Val column in CSV."
    if fdr_col != "FDR":
        df = df.rename(columns={fdr_col: "FDR"})

    # Optional PValue
    p_col = _pick_column(df, PVALUE_CANDIDATE_COLUMNS)
    if p_col and p_col != "PValue":
        df = df.rename(columns={p_col: "PValue"})

    # Optional logCPM / AveExpr
    logcpm_col = _pick_column(df, LOGCPM_CANDIDATE_COLUMNS)
    if logcpm_col and logcpm_col != "logCPM":
        df = df.rename(columns={logcpm_col: "logCPM"})

    ave_expr_col = _pick_column(df, AVE_EXPR_CANDIDATE_COLUMNS)
    if ave_expr_col and ave_expr_col != "AveExpr":
        df = df.rename(columns={ave_expr_col: "AveExpr"})

    # Standardize numeric types
    for col in ["logFC", "FDR", "PValue", "logCPM", "AveExpr"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Handle duplicates
    df = _handle_duplicates(df, duplicate_strategy)

    # Derived columns
    df["neglog10FDR"] = -np.log10(np.clip(df["FDR"].astype(float), 1e-300, None))
    if "PValue" in df.columns:
        df["neglog10PValue"] = -np.log10(np.clip(df["PValue"].astype(float), 1e-300, None))

    return df, None


@st.cache_data(show_spinner=False)
def load_all_contrasts_for_dataset(
    dataset_key: str,
    duplicate_strategy: str = "min_fdr",
) -> Dict[str, Tuple[Optional[pd.DataFrame], Optional[str]]]:
    """Pre-load all contrasts for a dataset. Returns mapping contrast → (df, error)."""
    ds_cfg = DATASETS.get(dataset_key, {})
    contrasts = ds_cfg.get("contrasts", {})
    result: Dict[str, Tuple[Optional[pd.DataFrame], Optional[str]]] = {}
    for contrast in contrasts:
        result[contrast] = load_and_standardize(dataset_key, contrast, duplicate_strategy=duplicate_strategy)
    return result


def get_available_contrasts(dataset_key: str) -> list[str]:
    ds_cfg = DATASETS.get(dataset_key, {})
    return list(ds_cfg.get("contrasts", {}).keys())

