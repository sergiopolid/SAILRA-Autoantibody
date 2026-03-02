## edgeR Autoantibody Explorer (Streamlit)

This app is a production-ready Streamlit dashboard for exploring edgeR differential autoantibody results across four datasets (HuSIGHT full length IgG/IgA and VirSIGHT IgG/IgA) plus a cross-isotype comparison view. It is designed to run locally and on Render.com with configuration- and environment-based data loading (no file uploads).

### Features

- **Four dataset tabs**
  - `HuSIGHT full length IgG`
  - `HuSIGHT full length IgA`
  - `VirSIGHT IgG`
  - `VirSIGHT IgA`
- **Per-dataset functionality**
  - Contrast selector (if multiple contrasts exist in `config.py`)
  - Interactive volcano plot (logFC vs \(-\log_{10}(\text{FDR})\)) with significance coloring and threshold lines
  - Interactive MA plot (logCPM or AveExpr vs logFC) with 0-line and friendly message if unavailable
  - Hits table with filtering (FDR, |logFC|, optional min logCPM, significant-only) and CSV export
  - Antigen search box + dropdown selection
  - Selected antigen details panel and highlighting across plots
- **Cross-isotype comparison tab**
  - `Compare IgG vs IgA` top-level tab
  - Select dataset family (`HuSIGHT full length` or `VirSIGHT`) and contrast (intersection of IgG/IgA contrasts)
  - IgG vs IgA logFC scatter with quadrant and cutoff lines
  - Side-by-side volcano and MA plots for IgG and IgA
  - Merged IgG/IgA table with `delta_logFC`, significance flags, and CSV export
  - Antigen selection that highlights across all comparison plots and the merged table

All interactions are scoped to the active tab/dataset, with per-dataset state stored in `st.session_state`.

---

### Configuration (`config.py`)

Datasets and contrasts are defined in `config.py`:

```python
DATASETS = {
    "husight_full_length_igg": {
        "label": "HuSIGHT full length IgG",
        "family": "husight_full_length",
        "isotype": "IgG",
        "contrasts": {
            "GroupA_vs_GroupB": {
                "csv": "./data/husight_full_length_igg_GroupA_vs_GroupB.csv",
            },
            # add more contrasts as needed
        },
    },
    # husight_full_length_iga, virsight_igg, virsight_iga ...
}
```

Each contrast supports either:

- **Local path**: `{"csv": "./data/xxx.csv"}`
- **Remote URL**: `{"csv_url": "https://.../xxx.csv"}`

`family` and `isotype` are used to wire up the `Compare IgG vs IgA` tab (families that have both IgG and IgA datasets are available for comparison).

---

### Environment variable overrides

CSV locations can be overridden per dataset (and optionally per contrast) using environment variables. This is Render-friendly and lets you avoid hard-coding file paths or URLs.

#### Naming scheme

The base for each dataset key is its uppercased, non-alphanumeric-normalized form. For example:

- Dataset key `husight_full_length_igg` → base `HUSIGHT_FULL_LENGTH_IGG`
- Dataset key `virsight_iga` → base `VIRSIGHT_IGA`

Contrast names are similarly normalized, e.g.:

- Contrast `GroupA_vs_GroupB` → `GROUPA_VS_GROUPB`

For each dataset and optional contrast, the app checks environment variables in this order:

1. `BASE_CONTRAST_CSV_URL`
2. `BASE_CONTRAST_CSV_PATH`
3. `BASE_CSV_URL`
4. `BASE_CSV_PATH`

Where:

- `BASE` is the normalized dataset key
- `CONTRAST` is the normalized contrast name

#### Examples

For dataset key `husight_full_length_igg` and contrast `GroupA_vs_GroupB`:

- `HUSIGHT_FULL_LENGTH_IGG_GROUPA_VS_GROUPB_CSV_URL`
- `HUSIGHT_FULL_LENGTH_IGG_GROUPA_VS_GROUPB_CSV_PATH`
- `HUSIGHT_FULL_LENGTH_IGG_CSV_URL`
- `HUSIGHT_FULL_LENGTH_IGG_CSV_PATH`

For dataset key `virsight_iga` and contrast `GroupA_vs_GroupB`:

- `VIRSIGHT_IGA_GROUPA_VS_GROUPB_CSV_URL`
- `VIRSIGHT_IGA_GROUPA_VS_GROUPB_CSV_PATH`
- `VIRSIGHT_IGA_CSV_URL`
- `VIRSIGHT_IGA_CSV_PATH`

**Behavior:**

- If an env var with a contrast suffix exists, it overrides the path/URL for that specific contrast.
- If only the base env var (without contrast suffix) exists, it serves as the **default contrast** location when that contrast is selected and no contrast-specific override is present.
- If neither `config.py` nor any env var resolves to a CSV, the UI shows a clear error message with the attempted env var names.

---

### CSV standardization and duplicate handling

Each CSV is expected to be an edgeR-style results table. The loader (`data.py`) robustly detects columns using the following preferences:

- **Antigen ID**: any of `["antigen","Antigen","gene","Gene","target","Target","symbol","Symbol","ID","id","Protein","protein","Name","name"]`
- **logFC (required)**: e.g. `"logFC"`, `"logfc"`, `"LogFC"`
- **FDR (required)**: any of `["FDR","fdr","adj.P.Val","adj.P.Val.","padj","qvalue","q_value"]`
- **PValue (optional)**: e.g. `"PValue"`, `"P.Value"`, `"p_value"`
- **logCPM (optional)**: any of `["logCPM","logcpm","LogCPM"]`
- **AveExpr (optional)**: any of `["AveExpr","ave_expr","Amean","avg_log_expr"]`

Columns are standardized to:

- `antigen`, `logFC`, `FDR`, and optionally `PValue`, `logCPM`, `AveExpr`
- Derived `neglog10FDR = -log10(max(FDR, 1e-300))`

Duplicate antigen handling strategies (user-selectable in the UI, per dataset or comparison):

- **Lowest FDR per antigen** (default; `min_fdr`)
- **Highest |logFC| per antigen** (`max_abs_logfc`)
- **Mean of numeric columns** (`mean`)

---

### UI overview

- **Per-dataset tabs** (`app.py` → `_render_dataset_tab`)
  - Contrast selector (if more than one contrast is configured)
  - FDR cutoff slider (`0–0.2`, default `0.05`)
  - |logFC| cutoff slider (`0–3`, default `0.7`)
  - Optional min logCPM cutoff slider (if `logCPM` present)
  - Toggle: show only significant hits (based on FDR + |logFC| thresholds)
  - Integer input: label top N points in the volcano
  - Ranking metric selector for labeling (`FDR`, `abs_logFC`, `neglog10FDR`)
  - Antigen search box and dropdown selection
  - Selected antigen details table (filtered to current view)
  - Interactive volcano and MA plots with selected-antigen highlighting
  - Hits table and a **download button** for the filtered CSV

- **Compare IgG vs IgA tab** (`app.py` → `_render_compare_tab`)
  - Family selector: `HuSIGHT full length` or `VirSIGHT` (families defined via `family` in `config.py`)
  - Contrast selector: intersection of IgG and IgA contrast names for the family
  - Shared FDR and |logFC| cutoffs
  - IgG vs IgA logFC scatter with:
    - Quadrant lines at 0
    - Optional cutoff lines at ±|logFC|
    - Categories: `Both`, `IgG only`, `IgA only`, `NS`
  - Side-by-side volcano plots for IgG and IgA
  - Side-by-side MA plots (if logCPM/AveExpr available; otherwise informative message)
  - Merged IgG/IgA table with:
    - `delta_logFC = logFC_IgA - logFC_IgG`
    - `sig_IgG`, `sig_IgA` boolean flags
  - Antigen search and selection that:
    - Highlights the selected antigen in scatter, volcano, and MA plots
    - Shows the joined IgG/IgA details in a dedicated panel
  - Download button for the merged CSV table

State is managed in `st.session_state` with keys that encode dataset, contrast, and family, so each tab’s interactions are scoped appropriately.

---

### Syncing with GitHub and Render

The repo is set up so **only the webapp** is tracked; large data and local reports are ignored.

**What gets committed (see `.gitignore`):**

- `app.py`, `config.py`, `data.py`, `plots.py`
- `requirements.txt`, `README.md`, `.gitignore`, `render.yaml`

**What is ignored:**

- `InfinityBio/` (all reports, scripts, and analysis outputs)
- `data/` and any `*.csv` (use environment variables or URLs for data on Render)
- Virtual envs, `__pycache__`, IDE files, `.env`, etc.

**Initial push to GitHub:**

```bash
cd /path/to/InfinityBio
git init
git add app.py config.py data.py plots.py requirements.txt README.md .gitignore render.yaml
git status   # confirm only these files are staged
git commit -m "Streamlit edgeR autoantibody dashboard for Render"
git branch -M main
git remote add origin https://github.com/YOUR_ORG/InfinityBio.git
git push -u origin main
```

After the first push, Render.com can connect to this GitHub repo. Configure the build/start commands and environment variables in the Render dashboard; no large files will be in the repo.

---

### Running locally

1. **Create and activate a virtual environment** (optional but recommended).
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Make sure your `config.py` paths/URLs (or environment variables) point to valid CSVs.
4. Launch Streamlit:

```bash
streamlit run app.py
```

5. Open the URL printed by Streamlit (usually `http://localhost:8501`) in a browser.

---

### Deploying on Render.com

1. **Create a new Web Service** in Render pointing at this repo.
2. Set **Environment** to Python.
3. Configure:
   - **Build command**:

     ```bash
     pip install -r requirements.txt
     ```

   - **Start command**:

     ```bash
     streamlit run app.py --server.port=$PORT --server.address=0.0.0.0
     ```

4. In the Render **Environment** settings, define any CSV override variables you need, e.g.:

   - `HUSIGHT_FULL_LENGTH_IGG_CSV_URL=https://.../husight_full_length_igg.csv`
   - `HUSIGHT_FULL_LENGTH_IGA_GROUPA_VS_GROUPB_CSV_URL=https://.../husight_full_length_iga_GroupA_vs_GroupB.csv`
   - `VIRSIGHT_IGG_CSV_PATH=/data/virsight_igg.csv`

5. Deploy; Render will build and run the app. Once live, you can adjust environment variables to repoint datasets without changing code.

---

### Notes and assumptions

- The app does **not** provide a file upload UI; all CSVs must be accessible via paths or URLs defined in `config.py` or environment variables.
- The MA plot is shown only when either `logCPM` or `AveExpr` is present (case/alias-insensitive); otherwise, a friendly “MA plot unavailable” message is displayed.
- Click-based selection from plots is not used due to current Streamlit limitations; selection is handled via search and dropdowns while still providing consistent cross-plot highlighting for the chosen antigen.

