#!/usr/bin/env python3
"""Export the notebook's PubMed search categories and construction commands."""

from __future__ import annotations

import ast
import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "pubmed_neuro_oncology_neuroscience_trends.ipynb"
OUTPUT_DIR = ROOT / "search_term_exports"


def collect_notebook_definitions():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    definitions = {}
    for cell_index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for node in tree.body:
            names = []
            if isinstance(node, ast.Assign):
                names = [
                    target.id for target in node.targets
                    if isinstance(target, ast.Name)
                ]
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names = [node.target.id]
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                names = [node.name]
            for name in names:
                definitions[name] = {
                    "source": ast.get_source_segment(source, node) or ast.unparse(node),
                    "node": node,
                    "cell": cell_index,
                }
    return definitions


def literal_string(definitions, name):
    node = definitions[name]["node"]
    value = node.value
    if isinstance(value, ast.Constant) and isinstance(value.value, str):
        return value.value
    if (
        isinstance(value, ast.Call)
        and isinstance(value.func, ast.Attribute)
        and value.func.attr == "strip"
        and isinstance(value.func.value, ast.Constant)
        and isinstance(value.func.value.value, str)
    ):
        return value.func.value.value.strip()
    raise ValueError(f"{name} is not a literal string assignment")


def query_section(definitions, name, heading=None):
    title = heading or name
    return (
        f"{title}\n"
        f"{'-' * len(title)}\n"
        f"{literal_string(definitions, name)}\n"
    )


def source_section(definitions, name, heading=None):
    title = heading or name
    definition = definitions[name]
    return (
        f"{title}\n"
        f"{'-' * len(title)}\n"
        f"Notebook cell: {definition['cell']}\n\n"
        f"{definition['source'].rstrip()}\n"
    )


def document(title, sections, notes=None):
    header = (
        f"{title}\n"
        f"{'=' * len(title)}\n\n"
        f"Source notebook: {NOTEBOOK.name}\n"
        "Publication interval used by the notebook: 2000-2025\n"
    )
    if notes:
        header += f"\nNotes\n-----\n{notes.strip()}\n"
    return header + "\n\n".join(section.rstrip() for section in sections) + "\n"


def main():
    definitions = collect_notebook_definitions()
    OUTPUT_DIR.mkdir(exist_ok=True)

    context_commands = [
        source_section(definitions, "NEURO_ONCOLOGY_DISEASE_CONTEXT"),
        source_section(definitions, "PRIMARY_CNS_DISEASE_CONTEXT"),
        source_section(definitions, "CNS_METASTASIS_DISEASE_CONTEXT"),
    ]

    molecular_target_sections = [
        source_section(definitions, "MOLECULAR_TARGET_FAMILIES"),
        query_section(definitions, "PRIMARY_LITERATURE_FILTER"),
        source_section(definitions, "build_target_query"),
        source_section(definitions, "get_target_year_counts"),
        source_section(definitions, "get_molecular_target_trends"),
    ]

    files = [
        (
            "01_neuro_oncology_terms.txt",
            document(
                "1. Neuro-Oncology Terms",
                [query_section(definitions, "NEURO_ONCOLOGY_QUERY")],
                "This is the broad neuro-oncology disease query used for the overall publication trend.",
            ),
        ),
        (
            "02_general_neurology_terms.txt",
            document(
                "2. General Neurology Terms",
                [query_section(definitions, "GENERAL_NEUROSCIENCE_QUERY")],
                "The notebook variable is named GENERAL_NEUROSCIENCE_QUERY.",
            ),
        ),
        (
            "03_small_molecule_and_targeted_therapy.txt",
            document(
                "3. Small-Molecule and Targeted Therapy",
                [
                    query_section(definitions, "SMALL_MOLECULE_TARGETED_THERAPY_QUERY"),
                    source_section(definitions, "SMALL_MOLECULE_NEURO_ONCOLOGY_QUERY"),
                    source_section(definitions, "NEURO_ONCOLOGY_DISEASE_CONTEXT"),
                ],
                "The first section is the therapy-term block; the following commands combine it with the neuro-oncology context.",
            ),
        ),
        (
            "04_antibody_and_adc.txt",
            document(
                "4. Antibody and ADC",
                [
                    query_section(definitions, "BIOLOGIC_ANTIBODY_THERAPY_QUERY"),
                    source_section(definitions, "BIOLOGIC_NEURO_ONCOLOGY_QUERY"),
                    source_section(definitions, "NEURO_ONCOLOGY_DISEASE_CONTEXT"),
                ],
            ),
        ),
        (
            "05_immunotherapy.txt",
            document(
                "5. Immunotherapy",
                [
                    query_section(definitions, "IMMUNOTHERAPY_QUERY"),
                    source_section(definitions, "IMMUNOTHERAPY_NEURO_ONCOLOGY_QUERY"),
                    source_section(definitions, "NEURO_ONCOLOGY_DISEASE_CONTEXT"),
                ],
            ),
        ),
        (
            "06_phase_iii_criteria.txt",
            document(
                "6. Phase III Criteria",
                [
                    query_section(definitions, "PHASE_III_CRITERIA"),
                    query_section(definitions, "PHASE_III_ARTICLE_TYPE_CRITERIA"),
                    query_section(definitions, "TREATMENT_PUBLICATION_EXCLUSION"),
                    query_section(definitions, "PHASE_III_CNS_METASTASIS_CONTEXT"),
                    source_section(definitions, "TREATMENT_MODALITY_QUERIES"),
                    *context_commands,
                    source_section(definitions, "get_phase_iii_context_modality_trends"),
                ],
                "Includes the clinical-trial article restriction, review/meta-analysis exclusions, radiation exclusion, analysis contexts, and final query builder.",
            ),
        ),
        (
            "07_molecular_targets_small_molecules.txt",
            document(
                "7. Molecular Targets for Small Molecules",
                molecular_target_sections,
                "Includes target families, aliases, primary-literature filtering, and query/count construction.",
            ),
        ),
        (
            "08_molecular_targets_small_molecules_2.txt",
            document(
                "8. Molecular Targets for Small Molecules",
                molecular_target_sections,
                "This category was listed twice in the request; this separately numbered export intentionally repeats item 7.",
            ),
        ),
        (
            "09_molecular_targets_biologics.txt",
            document(
                "9. Molecular Targets of Biologics",
                [
                    source_section(definitions, "BIOLOGIC_TARGET_FAMILIES"),
                    query_section(definitions, "PRIMARY_LITERATURE_FILTER"),
                    source_section(definitions, "build_target_query"),
                    source_section(definitions, "get_target_year_counts"),
                    source_section(definitions, "get_biologic_target_trends"),
                ],
                "Includes biologic target families, aliases, primary-literature filtering, and query/count construction.",
            ),
        ),
        (
            "10_molecular_targets_immunotherapy.txt",
            document(
                "10. Molecular Targets of Immunotherapy",
                [
                    source_section(definitions, "IMMUNE_CELL_THERAPY_FAMILIES"),
                    source_section(definitions, "IMMUNE_CELL_THERAPY_QUALIFIERS"),
                    query_section(definitions, "PRIMARY_LITERATURE_FILTER"),
                    source_section(definitions, "build_target_query"),
                    source_section(definitions, "build_immune_cell_therapy_query"),
                    source_section(definitions, "get_therapy_query_year_counts"),
                    source_section(definitions, "get_immune_cell_therapy_trends"),
                ],
                "Includes checkpoint, vaccine, CAR-T, and oncolytic-virus targets; therapy qualifiers; aliases; and query/count construction.",
            ),
        ),
    ]

    index_lines = [
        "SEARCH TERM EXPORT INDEX",
        "========================",
        "",
        f"Source notebook: {NOTEBOOK.name}",
        "",
    ]
    for filename, contents in files:
        (OUTPUT_DIR / filename).write_text(contents, encoding="utf-8")
        first_line = contents.splitlines()[0]
        index_lines.append(f"{filename}\n  {first_line}")
    index_lines.extend(
        [
            "",
            "NOTE",
            "----",
            "Items 7 and 8 use the same category wording in the supplied request.",
            "Both separately numbered files are included; item 8 intentionally repeats item 7.",
            "",
        ]
    )
    (OUTPUT_DIR / "00_index.txt").write_text(
        "\n".join(index_lines), encoding="utf-8"
    )

    archive = ROOT / "search_term_exports"
    shutil.make_archive(str(archive), "zip", ROOT, OUTPUT_DIR.name)
    print(f"Exported {len(files)} category TXT files to {OUTPUT_DIR}")
    print(f"Created archive: {archive.with_suffix('.zip')}")


if __name__ == "__main__":
    main()
