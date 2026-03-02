from __future__ import annotations

"""
Configuration for edgeR differential autoantibody Streamlit dashboard.

Datasets and contrasts can be customized here. Each contrast supports either:
- A local CSV path via the "csv" key
- A remote CSV URL via the "csv_url" key

Environment variables can override these locations. See README for details.
"""

DATASETS = {
    "husight_full_length_igg": {
        "label": "HuSIGHT full length IgG",
        "family": "husight_full_length",
        "isotype": "IgG",
        "contrasts": {
            "GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA": {
                # Copy the source CSV into ./webdata/ with this name:
                # /n/data1/bwh/medicine/kim/lab/lung/sergio/2026/InfinityBio/InfinityBio/BWH_SPoli_IB1251_IgG_Cohort_1_HuSIGHT_FullLength_IgG_Reports/edger_Ig_offset_results/edger_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA_Ig_offset.csv
                "csv": "./webdata/edger_HuSIGHT_FullLength_IgG_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA_Ig_offset.csv",
            },
        },
    },
    "husight_full_length_iga": {
        "label": "HuSIGHT full length IgA",
        "family": "husight_full_length",
        "isotype": "IgA",
        "contrasts": {
            "GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA": {
                # Copy the source CSV into ./webdata/ with this name:
                # /n/data1/bwh/medicine/kim/lab/lung/sergio/2026/InfinityBio/InfinityBio/BWH_SPoli_IB1251_IgA_Cohort_1_HuSIGHT_FullLength_IgA_Reports/edger_Ig_offset_results/edger_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA_Ig_offset.csv
                "csv": "./webdata/edger_HuSIGHT_FullLength_IgA_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA_Ig_offset.csv",
            },
        },
    },
    "virsight_igg": {
        "label": "VirSIGHT IgG",
        "family": "virsight",
        "isotype": "IgG",
        "contrasts": {
            "GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA": {
                # Copy the source CSV into ./webdata/ with this name:
                # /n/data1/bwh/medicine/kim/lab/lung/sergio/2026/InfinityBio/InfinityBio/BWH_SPoli_IB1251_IgG_Cohort_1_VirSIGHT_IgG_Reports/edger_Ig_offset_results/edger_VirSIGHT_IgG_species_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA.csv
                "csv": "./webdata/edger_VirSIGHT_IgG_species_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA.csv",
            },
        },
    },
    "virsight_iga": {
        "label": "VirSIGHT IgA",
        "family": "virsight",
        "isotype": "IgA",
        "contrasts": {
            "GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA": {
                # Copy the source CSV into ./webdata/ with this name:
                # /n/data1/bwh/medicine/kim/lab/lung/sergio/2026/InfinityBio/InfinityBio/BWH_SPoli_IB1251_IgA_Cohort_1_VirSIGHT_IgA_Reports/edger_Ig_offset_results/edger_VirSIGHT_IgA_species_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA.csv
                "csv": "./webdata/edger_VirSIGHT_IgA_species_GroupA_RA-ILD_RA-ILA_vs_GroupB_RA-noILA.csv",
            },
        },
    },
}


def get_dataset_keys_by_family() -> dict[str, dict[str, str]]:
    """
    Helper mapping of dataset family → {isotype → dataset_key}.

    Example:
    {
        "husight_full_length": {"IgG": "husight_full_length_igg", "IgA": "husight_full_length_iga"},
        "virsight": {"IgG": "virsight_igg", "IgA": "virsight_iga"},
    }
    """
    families: dict[str, dict[str, str]] = {}
    for key, cfg in DATASETS.items():
        family = cfg.get("family")
        isotype = cfg.get("isotype")
        if not family or not isotype:
            continue
        families.setdefault(family, {})[isotype] = key
    return families

