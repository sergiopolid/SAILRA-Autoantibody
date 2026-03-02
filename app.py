from __future__ import annotations

import io
from typing import Optional

import numpy as np
import pandas as pd
import streamlit as st

import config
from config import DATASETS, get_dataset_keys_by_family
from data import (
    get_available_contrasts,
    load_and_standardize,
)
from plots import igg_iga_scatter, ma_plot, volcano_plot


st.set_page_config(page_title="edgeR Autoantibody Explorer", layout="wide")


def _get_selected_antigen_key(dataset_key: str, contrast: str) -> str:
    return f"selected_antigen__{dataset_key}__{contrast}"


def _get_selected_antigen_compare_key(family_key: str) -> str:
    return f"selected_antigen_compare__{family_key}"


def _download_button(label: str, df: pd.DataFrame, filename: str, key: str) -> None:
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
        key=key,
    )


NO_SELECTION_LABEL = "— No selection —"


def _antigen_search_and_select(
    df: pd.DataFrame,
    state_key: str,
    title: str = "Search antigen",
) -> Optional[str]:
    """
    Search box + dropdown-based antigen selection.
    First option is "No selection" so the plot displays without highlighting any protein.
    """
    if state_key not in st.session_state:
        st.session_state[state_key] = None

    search_term = st.text_input(title, key=f"{state_key}__search")
    antigens = df["antigen"].astype(str)
    if search_term:
        mask = antigens.str.contains(search_term, case=False, na=False)
        filtered = sorted(antigens[mask].unique().tolist())
    else:
        filtered = sorted(antigens.unique().tolist())

    if not filtered:
        st.info("No antigens match the current search/filter criteria.")
        st.selectbox(
            "Select antigen to display",
            [NO_SELECTION_LABEL],
            index=0,
            key=f"{state_key}__selectbox",
        )
        st.session_state[state_key] = None
        return None

    options = [NO_SELECTION_LABEL] + filtered
    default_index = 0  # No selection
    current = st.session_state[state_key]
    if current and current in filtered:
        default_index = filtered.index(current) + 1

    selected = st.selectbox(
        "Select antigen to display",
        options,
        index=default_index,
        key=f"{state_key}__selectbox",
    )
    value = None if selected == NO_SELECTION_LABEL else selected
    st.session_state[state_key] = value
    return value


def _antigen_details_panel(df: pd.DataFrame, selected_antigen: Optional[str]) -> None:
    st.subheader("Selected antigen details")
    if not selected_antigen:
        st.info("Select an antigen using the search box to see details.")
        return
    rows = df[df["antigen"] == selected_antigen]
    if rows.empty:
        st.warning("Selected antigen not present in current filtered data.")
        return
    st.dataframe(rows.reset_index(drop=True))


def _filtered_df(
    df: pd.DataFrame,
    fdr_cutoff: float,
    logfc_cutoff: float,
    min_logcpm: Optional[float],
    show_only_sig: bool,
) -> pd.DataFrame:
    filtered = df.copy()
    if min_logcpm is not None and "logCPM" in filtered.columns:
        filtered = filtered[filtered["logCPM"] >= min_logcpm]
    sig_mask = (filtered["FDR"] <= fdr_cutoff) & (filtered["logFC"].abs() >= logfc_cutoff)
    if show_only_sig:
        filtered = filtered[sig_mask]
    return filtered


def _render_dataset_tab(dataset_key: str) -> None:
    ds_cfg = DATASETS[dataset_key]
    label = ds_cfg.get("label", dataset_key)
    contrasts = get_available_contrasts(dataset_key)

    # Controls within tab (scoped by dataset + contrast)
    st.markdown(f"### {label}")

    if not contrasts:
        st.error("No contrasts configured for this dataset.")
        return

    if len(contrasts) == 1:
        contrast = contrasts[0]
    else:
        contrast = st.selectbox("Contrast", contrasts, key=f"{dataset_key}__contrast")

    # Duplicate handling strategy (per dataset, global across contrasts)
    dup_strategy = st.selectbox(
        "Duplicate antigen handling",
        options=[
            ("min_fdr", "Lowest FDR per antigen"),
            ("max_abs_logfc", "Highest |logFC| per antigen"),
            ("mean", "Mean of numeric columns"),
        ],
        format_func=lambda x: x[1],
        key=f"{dataset_key}__dup_strategy",
    )[0]

    df, err = load_and_standardize(dataset_key, contrast, duplicate_strategy=dup_strategy)
    if err:
        st.error(err)
        return

    # Sliders and toggles
    col1, col2, col3 = st.columns(3)
    with col1:
        fdr_cutoff = st.slider(
            "FDR cutoff",
            min_value=0.0,
            max_value=0.2,
            value=0.05,
            step=0.005,
            key=f"{dataset_key}__{contrast}__fdr",
        )
    with col2:
        logfc_cutoff = st.slider(
            "|logFC| cutoff",
            min_value=0.0,
            max_value=3.0,
            value=0.7,
            step=0.05,
            key=f"{dataset_key}__{contrast}__logfc",
        )
    with col3:
        min_logcpm_val = None
        if "logCPM" in df.columns:
            min_logcpm_val = st.slider(
                "Min logCPM cutoff (optional)",
                min_value=float(np.nanmin(df["logCPM"])),
                max_value=float(np.nanmax(df["logCPM"])),
                value=float(np.nanmin(df["logCPM"])),
                step=0.5,
                key=f"{dataset_key}__{contrast}__logcpm",
            )

    col4, col5, col6 = st.columns(3)
    with col4:
        show_only_sig = st.checkbox(
            "Show only significant",
            value=False,
            key=f"{dataset_key}__{contrast}__sig_only",
        )
    with col5:
        label_top_n = st.number_input(
            "Label top N points",
            min_value=0,
            max_value=200,
            value=0,
            step=1,
            key=f"{dataset_key}__{contrast}__topn",
        )
    with col6:
        ranking_metric = st.selectbox(
            "Ranking metric for labels",
            options=["FDR", "abs_logFC", "neglog10FDR"],
            key=f"{dataset_key}__{contrast}__ranking",
        )

    # Filtered data for current view
    filtered = _filtered_df(df, fdr_cutoff, logfc_cutoff, min_logcpm_val, show_only_sig)

    # Antigen search & selection
    state_key = _get_selected_antigen_key(dataset_key, contrast)
    with st.expander("Antigen selection and details", expanded=True):
        selected_antigen = _antigen_search_and_select(filtered, state_key)
        _antigen_details_panel(filtered, selected_antigen)

    # Plots
    st.subheader("Plots")
    col_left, col_right = st.columns(2)

    with col_left:
        st.markdown("**Volcano plot**")
        fig_volcano = volcano_plot(
            filtered,
            fdr_cutoff=fdr_cutoff,
            logfc_cutoff=logfc_cutoff,
            label_top_n=label_top_n,
            ranking_metric=ranking_metric,
            selected_antigen=selected_antigen,
        )
        st.plotly_chart(fig_volcano, use_container_width=True)

    with col_right:
        st.markdown("**MA plot**")
        ma_x = None
        if "logCPM" in df.columns and not df["logCPM"].isna().all():
            ma_x = "logCPM"
        elif "AveExpr" in df.columns and not df["AveExpr"].isna().all():
            ma_x = "AveExpr"

        if ma_x is None:
            st.info("MA plot unavailable: logCPM/AveExpr column not found.")
        else:
            fig_ma = ma_plot(filtered, x_col=ma_x, selected_antigen=selected_antigen)
            st.plotly_chart(fig_ma, use_container_width=True)

    # Table + export
    st.subheader("Hits table")
    st.dataframe(filtered.reset_index(drop=True))
    _download_button(
        "Download filtered table as CSV",
        filtered,
        filename=f"{dataset_key}_{contrast}_filtered.csv",
        key=f"{dataset_key}__{contrast}__download",
    )


def _render_compare_tab() -> None:
    st.markdown("### Compare IgG vs IgA")

    families = get_dataset_keys_by_family()
    family_label_map = {
        "husight_full_length": "HuSIGHT full length",
        "virsight": "VirSIGHT",
    }
    family_options = [f for f in families.keys() if set(families[f].keys()) >= {"IgG", "IgA"}]
    if not family_options:
        st.error("No dataset families with both IgG and IgA configured.")
        return

    family_key = st.selectbox(
        "Dataset family",
        options=family_options,
        format_func=lambda x: family_label_map.get(x, x),
        key="compare__family",
    )

    igg_ds_key = families[family_key]["IgG"]
    iga_ds_key = families[family_key]["IgA"]

    igg_contrasts = set(get_available_contrasts(igg_ds_key))
    iga_contrasts = set(get_available_contrasts(iga_ds_key))
    common_contrasts = sorted(igg_contrasts & iga_contrasts)
    if not common_contrasts:
        st.error("No common contrasts found between IgG and IgA for this family.")
        return

    contrast = st.selectbox(
        "Contrast",
        options=common_contrasts,
        key=f"compare__{family_key}__contrast",
    )

    # Duplicate strategy shared across IgG/IgA
    dup_strategy = st.selectbox(
        "Duplicate antigen handling",
        options=[
            ("min_fdr", "Lowest FDR per antigen"),
            ("max_abs_logfc", "Highest |logFC| per antigen"),
            ("mean", "Mean of numeric columns"),
        ],
        format_func=lambda x: x[1],
        key=f"compare__{family_key}__dup_strategy",
    )[0]

    igg_df, igg_err = load_and_standardize(igg_ds_key, contrast, duplicate_strategy=dup_strategy)
    iga_df, iga_err = load_and_standardize(iga_ds_key, contrast, duplicate_strategy=dup_strategy)
    if igg_err:
        st.error(f"IgG: {igg_err}")
        return
    if iga_err:
        st.error(f"IgA: {iga_err}")
        return

    # Shared cutoffs
    col1, col2 = st.columns(2)
    with col1:
        fdr_cutoff = st.slider(
            "FDR cutoff",
            min_value=0.0,
            max_value=0.2,
            value=0.05,
            step=0.005,
            key=f"compare__{family_key}__{contrast}__fdr",
        )
    with col2:
        logfc_cutoff = st.slider(
            "|logFC| cutoff",
            min_value=0.0,
            max_value=3.0,
            value=0.7,
            step=0.05,
            key=f"compare__{family_key}__{contrast}__logfc",
        )

    # Merge IgG/IgA by antigen
    merge_cols = ["antigen", "logFC", "FDR", "PValue", "logCPM", "AveExpr"]
    igg_prefixed = igg_df[[c for c in merge_cols if c in igg_df.columns]].copy()
    iga_prefixed = iga_df[[c for c in merge_cols if c in iga_df.columns]].copy()
    igg_prefixed = igg_prefixed.rename(
        columns={c: f"{c}_IgG" for c in igg_prefixed.columns if c != "antigen"},
    )
    iga_prefixed = iga_prefixed.rename(
        columns={c: f"{c}_IgA" for c in iga_prefixed.columns if c != "antigen"},
    )
    merged = pd.merge(igg_prefixed, iga_prefixed, on="antigen", how="inner")
    merged["delta_logFC"] = merged["logFC_IgA"] - merged["logFC_IgG"]

    # Significance flags
    merged["sig_IgG"] = (merged["FDR_IgG"] <= fdr_cutoff) & (merged["logFC_IgG"].abs() >= logfc_cutoff)
    merged["sig_IgA"] = (merged["FDR_IgA"] <= fdr_cutoff) & (merged["logFC_IgA"].abs() >= logfc_cutoff)

    # Antigen selection for compare
    state_key = _get_selected_antigen_compare_key(family_key)
    with st.expander("Antigen selection and details", expanded=True):
        selected_antigen = _antigen_search_and_select(merged.rename(columns={"antigen": "antigen"}), state_key)
        st.subheader("Selected antigen (IgG & IgA)")
        if selected_antigen:
            rows = merged[merged["antigen"] == selected_antigen]
            if rows.empty:
                st.warning("Selected antigen not present in merged table.")
            else:
                st.dataframe(rows.reset_index(drop=True))
        else:
            st.info("Select an antigen using the search box to see IgG vs IgA details.")

    # Plots
    st.subheader("IgG vs IgA plots")
    col_scatter, col_volcano = st.columns(2)
    with col_scatter:
        st.markdown("**IgG vs IgA logFC scatter**")
        fig_scatter = igg_iga_scatter(
            merged,
            fdr_cutoff=fdr_cutoff,
            logfc_cutoff=logfc_cutoff,
            selected_antigen=selected_antigen,
        )
        st.plotly_chart(fig_scatter, use_container_width=True)

    with col_volcano:
        st.markdown("**Volcano plots (side-by-side)**")
        subcol1, subcol2 = st.columns(2)
        with subcol1:
            st.caption("IgG")
            fig_vol_igg = volcano_plot(
                igg_df,
                fdr_cutoff=fdr_cutoff,
                logfc_cutoff=logfc_cutoff,
                label_top_n=0,
                ranking_metric="FDR",
                selected_antigen=selected_antigen,
            )
            st.plotly_chart(fig_vol_igg, use_container_width=True)
        with subcol2:
            st.caption("IgA")
            fig_vol_iga = volcano_plot(
                iga_df,
                fdr_cutoff=fdr_cutoff,
                logfc_cutoff=logfc_cutoff,
                label_top_n=0,
                ranking_metric="FDR",
                selected_antigen=selected_antigen,
            )
            st.plotly_chart(fig_vol_iga, use_container_width=True)

    st.subheader("MA plots (side-by-side)")
    ma_col1, ma_col2 = st.columns(2)
    with ma_col1:
        st.caption("IgG")
        ma_x_igg = None
        if "logCPM" in igg_df.columns and not igg_df["logCPM"].isna().all():
            ma_x_igg = "logCPM"
        elif "AveExpr" in igg_df.columns and not igg_df["AveExpr"].isna().all():
            ma_x_igg = "AveExpr"
        if ma_x_igg is None:
            st.info("MA plot unavailable for IgG: logCPM/AveExpr column not found.")
        else:
            fig_ma_igg = ma_plot(igg_df, x_col=ma_x_igg, selected_antigen=selected_antigen)
            st.plotly_chart(fig_ma_igg, use_container_width=True)

    with ma_col2:
        st.caption("IgA")
        ma_x_iga = None
        if "logCPM" in iga_df.columns and not iga_df["logCPM"].isna().all():
            ma_x_iga = "logCPM"
        elif "AveExpr" in iga_df.columns and not iga_df["AveExpr"].isna().all():
            ma_x_iga = "AveExpr"
        if ma_x_iga is None:
            st.info("MA plot unavailable for IgA: logCPM/AveExpr column not found.")
        else:
            fig_ma_iga = ma_plot(iga_df, x_col=ma_x_iga, selected_antigen=selected_antigen)
            st.plotly_chart(fig_ma_iga, use_container_width=True)

    # Merged table + export
    st.subheader("Merged IgG/IgA table")
    st.dataframe(merged.reset_index(drop=True))
    _download_button(
        "Download merged table as CSV",
        merged,
        filename=f"{family_key}_{contrast}_igg_vs_iga_merged.csv",
        key=f"compare__{family_key}__{contrast}__download",
    )


def main() -> None:
    st.title("edgeR Autoantibody Differential Analysis Explorer")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "HuSIGHT full length IgG",
            "HuSIGHT full length IgA",
            "VirSIGHT IgG",
            "VirSIGHT IgA",
            "Compare IgG vs IgA",
        ]
    )

    with tab1:
        _render_dataset_tab("husight_full_length_igg")
    with tab2:
        _render_dataset_tab("husight_full_length_iga")
    with tab3:
        _render_dataset_tab("virsight_igg")
    with tab4:
        _render_dataset_tab("virsight_iga")
    with tab5:
        _render_compare_tab()


if __name__ == "__main__":
    main()

