from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def _significance_category(
    df: pd.DataFrame,
    fdr_cutoff: float,
    logfc_cutoff: float,
) -> pd.Series:
    sig = (df["FDR"] <= fdr_cutoff) & (df["logFC"].abs() >= logfc_cutoff)
    up = sig & (df["logFC"] > 0)
    down = sig & (df["logFC"] < 0)
    cat = np.full(len(df), "NS", dtype=object)
    cat[up.to_numpy()] = "Up"
    cat[down.to_numpy()] = "Down"
    return pd.Series(cat, index=df.index, name="category")


def volcano_plot(
    df: pd.DataFrame,
    fdr_cutoff: float,
    logfc_cutoff: float,
    label_top_n: int = 0,
    ranking_metric: str = "FDR",
    selected_antigen: Optional[str] = None,
):
    data = df.copy()
    data["category"] = _significance_category(data, fdr_cutoff, logfc_cutoff)

    # Ranking metric for labeling
    if ranking_metric == "FDR":
        score = data["FDR"].fillna(1.0)
        order = score
        ascending = True
    elif ranking_metric == "abs_logFC":
        score = data["logFC"].abs()
        order = -score
        ascending = False
    else:  # neglog10FDR
        score = data["neglog10FDR"]
        order = -score
        ascending = False
    data["_rank_score"] = score

    text = None
    if label_top_n > 0:
        top_idx = order.sort_values(ascending=ascending).index[:label_top_n]
        text = data["antigen"].where(data.index.isin(top_idx))

    fig = px.scatter(
        data,
        x="logFC",
        y="neglog10FDR",
        color="category",
        hover_data={
            "antigen": True,
            "logFC": True,
            "FDR": True,
            "PValue": "PValue" in data.columns,
        },
        text=text,
        color_discrete_map={"Up": "#d62728", "Down": "#1f77b4", "NS": "#bbbbbb"},
    )

    # Threshold lines
    y_thresh = -np.log10(max(fdr_cutoff, 1e-300))
    fig.add_hline(y=y_thresh, line_dash="dot", line_color="black")
    fig.add_vline(x=logfc_cutoff, line_dash="dot", line_color="black")
    fig.add_vline(x=-logfc_cutoff, line_dash="dot", line_color="black")

    # Highlight selected antigen
    if selected_antigen and selected_antigen in data["antigen"].values:
        sel = data[data["antigen"] == selected_antigen]
        if not sel.empty:
            fig.add_trace(
                go.Scatter(
                    x=sel["logFC"],
                    y=sel["neglog10FDR"],
                    mode="markers+text",
                    text=sel["antigen"],
                    textposition="top center",
                    marker=dict(size=14, color="gold", line=dict(width=2, color="black")),
                    name="Selected antigen",
                    showlegend=True,
                )
            )

    fig.update_layout(
        xaxis_title="logFC",
        yaxis_title="-log10(FDR)",
        legend_title="Significance",
        template="plotly_white",
    )
    return fig


def ma_plot(
    df: pd.DataFrame,
    x_col: str,
    selected_antigen: Optional[str] = None,
):
    if x_col not in df.columns:
        raise ValueError(f"Column {x_col} not found in DataFrame.")

    data = df.copy()
    fig = px.scatter(
        data,
        x=x_col,
        y="logFC",
        hover_data={
            "antigen": True,
            "logFC": True,
            "FDR": True,
            x_col: True,
        },
        color_discrete_sequence=["#1f77b4"],
    )

    fig.add_hline(y=0.0, line_dash="dash", line_color="black")

    if selected_antigen and selected_antigen in data["antigen"].values:
        sel = data[data["antigen"] == selected_antigen]
        if not sel.empty:
            fig.add_trace(
                go.Scatter(
                    x=sel[x_col],
                    y=sel["logFC"],
                    mode="markers+text",
                    text=sel["antigen"],
                    textposition="top center",
                    marker=dict(size=14, color="gold", line=dict(width=2, color="black")),
                    name="Selected antigen",
                    showlegend=False,
                )
            )

    fig.update_layout(
        xaxis_title=x_col,
        yaxis_title="logFC",
        template="plotly_white",
    )
    return fig


def igg_iga_scatter(
    merged_df: pd.DataFrame,
    fdr_cutoff: float,
    logfc_cutoff: float,
    selected_antigen: Optional[str] = None,
):
    data = merged_df.copy()
    # Simple significance classification for both isotypes
    sig_igg = (data["FDR_IgG"] <= fdr_cutoff) & (data["logFC_IgG"].abs() >= logfc_cutoff)
    sig_iga = (data["FDR_IgA"] <= fdr_cutoff) & (data["logFC_IgA"].abs() >= logfc_cutoff)
    cat = np.full(len(data), "NS", dtype=object)
    cat[sig_igg & sig_iga] = "Both"
    cat[sig_igg & ~sig_iga] = "IgG only"
    cat[~sig_igg & sig_iga] = "IgA only"
    data["category"] = cat

    fig = px.scatter(
        data,
        x="logFC_IgG",
        y="logFC_IgA",
        color="category",
        hover_data={
            "antigen": True,
            "logFC_IgG": True,
            "FDR_IgG": True,
            "logFC_IgA": True,
            "FDR_IgA": True,
        },
        color_discrete_map={
            "Both": "#9467bd",
            "IgG only": "#1f77b4",
            "IgA only": "#d62728",
            "NS": "#bbbbbb",
        },
    )

    # Quadrant lines
    fig.add_hline(y=0.0, line_dash="dash", line_color="black")
    fig.add_vline(x=0.0, line_dash="dash", line_color="black")

    # Optional cutoffs
    fig.add_hline(y=logfc_cutoff, line_dash="dot", line_color="gray")
    fig.add_hline(y=-logfc_cutoff, line_dash="dot", line_color="gray")
    fig.add_vline(x=logfc_cutoff, line_dash="dot", line_color="gray")
    fig.add_vline(x=-logfc_cutoff, line_dash="dot", line_color="gray")

    if selected_antigen and selected_antigen in data["antigen"].values:
        sel = data[data["antigen"] == selected_antigen]
        if not sel.empty:
            fig.add_trace(
                go.Scatter(
                    x=sel["logFC_IgG"],
                    y=sel["logFC_IgA"],
                    mode="markers+text",
                    text=sel["antigen"],
                    textposition="top center",
                    marker=dict(size=14, color="gold", line=dict(width=2, color="black")),
                    name="Selected antigen",
                    showlegend=True,
                )
            )

    fig.update_layout(
        xaxis_title="logFC (IgG)",
        yaxis_title="logFC (IgA)",
        legend_title="Significance",
        template="plotly_white",
    )
    return fig

