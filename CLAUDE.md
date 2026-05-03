# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the App

All commands should be run from the `SDM-Primer-Designer-main/` directory (the project root):

```bash
# Install dependencies
pip install -r requirements.txt

# Run the Streamlit app
streamlit run main.py
```

The app runs on port 8501 by default. For Codespaces/devcontainer, CORS and XSRF protection are disabled:
```bash
streamlit run main.py --server.enableCORS false --server.enableXsrfProtection false
```

## Architecture

The codebase has two layers:

**`sdm_designer.py` — Core logic (no Streamlit dependency)**
- `SDMPrimerDesigner`: takes a template DNA sequence at init. Key methods:
  - `design(row, method, target_tm)`: designs forward + reverse primers for a mutation row (dict with `mutation_name`, `mode`, `aa_pos`, `target_aa`/`insert_seq`/`del_len`). Uses overlapping extension method, iterating primer arm lengths (15–44 nt) until Nearest-Neighbor Tm ≥ `target_tm`. Returns a result dict with `fwd_primer`, `rev_primer`, `Tm`, `GC_%`, `MW`, `Dimer_Tm`, `Product_Size(bp)`, `New_Sites`, `full_seq`, `mut_start`, `mut_end`.
  - `detect_features(sequence, custom_library)`: scans for known plasmid parts (AmpR, ori, T7/CMV promoters) plus any user-supplied custom parts.
  - `_calculate_self_dimer_tm(seq)`: sliding-window NN Tm check for self-complementarity.
  - `generate_protocol(result)`: wraps a result dict into a `PCRProtocol`.
- `PCRProtocol` (dataclass): computes PCR conditions at `__post_init__`: extension time (5 sec/kb), annealing time (Wallace rule vs. 55 °C threshold), water volume. Produces reaction tables, cycling text, and full downloadable text via `to_reaction_table()`, `to_cycling_text()`, `to_full_text()`.

**`main.py` — Streamlit UI**
- Sidebar: FASTA upload, mutation list (CSV/XLSX) upload, Tm slider, map view mode, custom feature registration/JSON export.
- Main flow (triggered by button): reads FASTA → `SeqIO.read` → `SDMPrimerDesigner` → loops `design()` per row → builds result DataFrame → renders styled table (Dimer_Tm heat coloring), Excel download (with embedded map images via `xlsxwriter`), and order-format text area.
- Post-run sections (if `st.session_state['results']` exists): primer dissolution guide, PCR protocol tabs (bulk download / per-mutation), vector map preview.

## Input File Formats

- FASTA: single-record `.fasta`/`.fa` (see `examples/template.fasta`)
- Mutation CSV columns: `mutation_name`, `mode` (`sub`/`ins`/`del`), `aa_pos`, `target_aa`, `insert_seq`, `del_len` (see `examples/mutations.csv`)

## Key Dependencies

| Package | Role |
|---|---|
| `biopython` | `Tm_NN`, `molecular_weight`, `Restriction`, `SeqIO` |
| `dna_features_viewer` | Linear/Circular vector map rendering |
| `xlsxwriter` | Excel report with embedded images |
| `streamlit` | UI framework |

## Notes

- The `tests/` directory exists but contains only `__init__.py`; no test runner is configured.
- Codon table (`ECOLI_CODONS`) and feature signatures (`FEATURE_LIBRARY`) are class-level constants on `SDMPrimerDesigner` and can be extended without breaking the interface.
