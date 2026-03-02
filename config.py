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
            # Example placeholder contrast; replace paths/URLs with real ones.
            "GroupA_vs_GroupB": {
                "csv": "./data/husight_full_length_igg_GroupA_vs_GroupB.csv",
            },
        },
    },
    "husight_full_length_iga": {
        "label": "HuSIGHT full length IgA",
        "family": "husight_full_length",
        "isotype": "IgA",
        "contrasts": {
            "GroupA_vs_GroupB": {
                "csv": "./data/husight_full_length_iga_GroupA_vs_GroupB.csv",
            },
        },
    },
    "virsight_igg": {
        "label": "VirSIGHT IgG",
        "family": "virsight",
        "isotype": "IgG",
        "contrasts": {
            "GroupA_vs_GroupB": {
                "csv": "./data/virsight_igg_GroupA_vs_GroupB.csv",
            },
        },
    },
    "virsight_iga": {
        "label": "VirSIGHT IgA",
        "family": "virsight",
        "isotype": "IgA",
        "contrasts": {
            "GroupA_vs_GroupB": {
                "csv": "./data/virsight_iga_GroupA_vs_GroupB.csv",
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

