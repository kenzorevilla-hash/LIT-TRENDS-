#!/usr/bin/env python3
"""Export the data plotted by every graph in the literature-trends notebook."""

from __future__ import annotations

import contextlib
import io
import json
import os
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "pubmed_neuro_oncology_neuroscience_trends.ipynb"
OUTPUT_DIR = ROOT / "graph_counts_csv"

# Keep matplotlib headless and prevent it from writing configuration outside the workspace.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "lit_trends_matplotlib")
)


def load_notebook_namespace() -> dict:
    """Run notebook code cells so exports use the notebook's exact cached datasets."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    namespace = {
        "__name__": "__notebook_export__",
        "display": lambda *args, **kwargs: None,
    }

    import matplotlib.pyplot as plt

    for cell_index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        if source.lstrip().startswith(("!", "%")):
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                exec(compile(source, f"{NOTEBOOK.name}:cell-{cell_index}", "exec"), namespace)
        except Exception as exc:
            first_line = source.strip().splitlines()[0]
            raise RuntimeError(
                f"Notebook cell {cell_index} failed while exporting ({first_line!r})"
            ) from exc
        finally:
            plt.close("all")

    return namespace


def integerize_counts(frame, columns):
    result = frame.copy()
    for column in columns:
        result[column] = result[column].astype("int64")
    return result


def cumulative_target_data(frame):
    result = frame.copy().sort_values(["Family", "Target", "Year"])
    result["Cumulative publications"] = result.groupby(
        ["Family", "Target"], sort=False
    )["Publications"].cumsum()
    result = integerize_counts(
        result,
        ["Year", "Publications", "Cumulative publications"],
    )
    return result[
        [
            "Context",
            "Family",
            "Target",
            "Year",
            "Publications",
            "Cumulative publications",
        ]
    ]


def main() -> None:
    namespace = load_notebook_namespace()
    OUTPUT_DIR.mkdir(exist_ok=True)

    df_pubs = namespace["df_pubs"].copy()
    df_pubs = integerize_counts(
        df_pubs,
        ["Year", "General neuroscience", "Neuro-oncology"],
    )

    share = df_pubs[["Year", "Neuro-oncology proportion"]].copy()
    share["Neuro-oncology share (%)"] = share["Neuro-oncology proportion"] * 100

    graph_exports = [
        (
            "01_neuro_oncology_publications_annual.csv",
            df_pubs[["Year", "Neuro-oncology"]],
            "Annual neuro-oncology publication counts",
            "Annual",
        ),
        (
            "02_general_neuroscience_vs_neuro_oncology_annual.csv",
            df_pubs[["Year", "General neuroscience", "Neuro-oncology"]],
            "Annual general-neuroscience and neuro-oncology publication counts",
            "Annual",
        ),
        (
            "03_neuro_oncology_share_annual.csv",
            share,
            "Annual neuro-oncology share of general-neuroscience publications",
            "Annual ratio and percent",
        ),
    ]

    treatment_exports = [
        (
            "04_neuro_oncology_treatment_modalities_cumulative.csv",
            "combined_treatment_cumulative",
            "Cumulative neuro-oncology treatment publications by modality",
        ),
        (
            "05_primary_cns_treatment_modalities_cumulative.csv",
            "primary_cns_treatment_cumulative",
            "Cumulative primary-CNS treatment publications by modality",
        ),
        (
            "06_cns_metastasis_treatment_modalities_cumulative.csv",
            "cns_metastasis_treatment_cumulative",
            "Cumulative CNS-metastasis treatment publications by modality",
        ),
        (
            "07_phase_iii_neuro_oncology_cumulative.csv",
            "phase_iii_overall_cumulative",
            "Cumulative Phase III neuro-oncology publications by modality",
        ),
        (
            "08_phase_iii_primary_cns_cumulative.csv",
            "phase_iii_primary_cns_cumulative",
            "Cumulative Phase III primary-CNS publications by modality",
        ),
        (
            "09_phase_iii_cns_metastasis_cumulative.csv",
            "phase_iii_cns_metastasis_cumulative",
            "Cumulative Phase III CNS-metastasis publications by modality",
        ),
    ]
    for filename, variable, description in treatment_exports:
        frame = namespace[variable].copy().sort_values(["Category", "Year"])
        frame = integerize_counts(
            frame,
            ["Year", "Publications", "Cumulative publications"],
        )
        graph_exports.append(
            (
                filename,
                frame[["Category", "Year", "Publications", "Cumulative publications"]],
                description,
                "Annual and cumulative",
            )
        )

    target_exports = [
        (
            "10_molecular_targets_neuro_oncology_cumulative.csv",
            "all_neuro_oncology_target_trends",
            "Cumulative molecular-target publications in neuro-oncology",
        ),
        (
            "11_molecular_targets_primary_cns_cumulative.csv",
            "primary_cns_target_trends",
            "Cumulative molecular-target publications in primary CNS tumors",
        ),
        (
            "12_molecular_targets_cns_metastasis_cumulative.csv",
            "cns_metastasis_target_trends",
            "Cumulative molecular-target publications in CNS metastasis",
        ),
        (
            "13_biologic_targets_neuro_oncology_cumulative.csv",
            "all_neuro_oncology_biologic_target_trends",
            "Cumulative biologic-target publications in neuro-oncology",
        ),
        (
            "14_biologic_targets_primary_cns_cumulative.csv",
            "primary_cns_biologic_target_trends",
            "Cumulative biologic-target publications in primary CNS tumors",
        ),
        (
            "15_biologic_targets_cns_metastasis_cumulative.csv",
            "cns_metastasis_biologic_target_trends",
            "Cumulative biologic-target publications in CNS metastasis",
        ),
        (
            "16_immune_cell_therapy_neuro_oncology_cumulative.csv",
            "all_neuro_oncology_immune_cell_trends",
            "Cumulative immune/cell-therapy publications in neuro-oncology",
        ),
        (
            "17_immune_cell_therapy_primary_cns_cumulative.csv",
            "primary_cns_immune_cell_trends",
            "Cumulative immune/cell-therapy publications in primary CNS tumors",
        ),
        (
            "18_immune_cell_therapy_cns_metastasis_cumulative.csv",
            "cns_metastasis_immune_cell_trends",
            "Cumulative immune/cell-therapy publications in CNS metastasis",
        ),
    ]
    for filename, variable, description in target_exports:
        graph_exports.append(
            (
                filename,
                cumulative_target_data(namespace[variable]),
                description,
                "Annual and cumulative",
            )
        )

    keyword_trends = namespace["keyword_trends"].copy().sort_values(["Rank", "Year"])
    keyword_trends = integerize_counts(
        keyword_trends,
        ["Rank", "Year", "Publications"],
    )
    graph_exports.append(
        (
            "19_top_10_2025_disease_keywords_annual_trends.csv",
            keyword_trends[["Rank", "Keyword", "Year", "Publications"]],
            "Annual trends for the ten leading 2025 prespecified neuro-oncology disease keywords",
            "Annual",
        )
    )

    index_rows = []
    for graph_number, (filename, frame, description, count_type) in enumerate(
        graph_exports, start=1
    ):
        frame.to_csv(OUTPUT_DIR / filename, index=False)
        index_rows.append(
            {
                "Graph number": graph_number,
                "CSV file": filename,
                "Description": description,
                "Count type": count_type,
                "Rows": len(frame),
            }
        )

    import pandas as pd

    pd.DataFrame(index_rows).to_csv(OUTPUT_DIR / "00_graph_counts_index.csv", index=False)
    print(f"Exported {len(graph_exports)} graph CSV files to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
