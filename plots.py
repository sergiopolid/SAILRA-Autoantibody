from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def _significance_category(
    df: pd.DataFrame,
    p_col: str,
    p_cutoff: float,
    logfc_cutoff: float,
) -> pd.Series:
    sig = (df[p_col] <= p_cutoff) & (df["logFC"].abs() >= logfc_cutoff)
    up = sig & (df["logFC"] > 0)
    down = sig & (df["logFC"] < 0)
    cat = np.full(len(df), "NS", dtype=object)
    cat[up.to_numpy()] = "Up"
    cat[down.to_numpy()] = "Down"
    return pd.Series(cat, index=df.index, name="category")


def volcano_plot(
    df: pd.DataFrame,
    p_col: str,
    p_cutoff: float,
    logfc_cutoff: float,
    label_top_n: int = 0,
    ranking_metric: str = "FDR",
    selected_antigen: Optional[str] = None,
):
    """p_col: 'FDR' or 'PValue'. Uses corresponding neglog10 column for y-axis."""
    data = df.copy()
    data["category"] = _significance_category(data, p_col, p_cutoff, logfc_cutoff)

    neglog10_col = "neglog10FDR" if p_col == "FDR" else "neglog10PValue"
    if neglog10_col not in data.columns:
        neglog10_col = "neglog10FDR"

    # Ranking metric for labeling
    if ranking_metric == "FDR":
        score = data["FDR"].fillna(1.0)
        order = score
        ascending = True
    elif ranking_metric == "PValue" and "PValue" in data.columns:
        score = data["PValue"].fillna(1.0)
        order = score
        ascending = True
    elif ranking_metric == "abs_logFC":
        score = data["logFC"].abs()
        order = -score
        ascending = False
    elif ranking_metric == "neglog10PValue" and "neglog10PValue" in data.columns:
        score = data["neglog10PValue"]
        order = -score
        ascending = False
    else:
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
        y=neglog10_col,
        color="category",
        hover_data={
            "antigen": True,
            "logFC": True,
            "FDR": "FDR" in data.columns,
            "PValue": "PValue" in data.columns,
        },
        text=text,
        color_discrete_map={"Up": "#d62728", "Down": "#1f77b4", "NS": "#bbbbbb"},
    )

    # Threshold lines
    y_thresh = -np.log10(max(p_cutoff, 1e-300))
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
                    y=sel[neglog10_col],
                    mode="markers+text",
                    text=sel["antigen"],
                    textposition="top center",
                    marker=dict(size=14, color="gold", line=dict(width=2, color="black")),
                    name="Selected antigen",
                    showlegend=True,
                )
            )

    y_label = "-log10(FDR)" if p_col == "FDR" else "-log10(P-value)"
    fig.update_layout(
        xaxis_title="logFC",
        yaxis_title=y_label,
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
    p_col_igg: str,
    p_col_iga: str,
    p_cutoff: float,
    logfc_cutoff: float,
    selected_antigen: Optional[str] = None,
):
    """p_col_igg/iga: e.g. 'FDR_IgG'/'FDR_IgA' or 'PValue_IgG'/'PValue_IgA'."""
    data = merged_df.copy()
    if p_col_igg not in data.columns:
        p_col_igg = "FDR_IgG"
    if p_col_iga not in data.columns:
        p_col_iga = "FDR_IgA"
    sig_igg = (data[p_col_igg] <= p_cutoff) & (data["logFC_IgG"].abs() >= logfc_cutoff)
    sig_iga = (data[p_col_iga] <= p_cutoff) & (data["logFC_IgA"].abs() >= logfc_cutoff)
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

