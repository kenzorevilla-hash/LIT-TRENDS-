#!/usr/bin/env python3
"""Export the PubMed studies underlying the phase III CNS-metastasis figure."""

from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import re
import sqlite3
import time
from urllib.error import HTTPError, URLError
from xml.sax.saxutils import escape, quoteattr
import zipfile

from Bio import Entrez


ROOT = Path(__file__).resolve().parent
NOTEBOOK_PATH = ROOT / "pubmed_neuro_oncology_neuroscience_trends.ipynb"
CACHE_PATH = ROOT / ".pubmed_cache" / "pubmed_counts.sqlite3"
DEFAULT_OUTPUT = ROOT / "phase_III_CNS_metastasis_studies_2000_2025.xlsx"
START_YEAR = 2000
END_YEAR = 2025


def _static_string(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if (
        isinstance(node, ast.Call)
        and not node.args
        and not node.keywords
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "strip"
    ):
        value = _static_string(node.func.value)
        return value.strip() if value is not None else None
    return None


def extract_notebook_queries() -> tuple[str, dict[str, str], str, str, str]:
    """Extract the exact query constants used by the plotted notebook cells."""
    with NOTEBOOK_PATH.open(encoding="utf-8") as handle:
        notebook = json.load(handle)

    strings: dict[str, str] = {}
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for statement in tree.body:
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
                continue
            target = statement.targets[0]
            if not isinstance(target, ast.Name):
                continue
            value = _static_string(statement.value)
            if value is not None:
                strings[target.id] = value

    required = {
        "CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY",
        "SMALL_MOLECULE_TARGETED_THERAPY_QUERY",
        "BIOLOGIC_ANTIBODY_THERAPY_QUERY",
        "IMMUNOTHERAPY_QUERY",
        "PHASE_III_CNS_METASTASIS_CONTEXT",
        "PHASE_III_CRITERIA",
        "PHASE_III_ARTICLE_TYPE_CRITERIA",
        "TREATMENT_PUBLICATION_EXCLUSION",
    }
    missing = sorted(required.difference(strings))
    if missing:
        raise RuntimeError(f"Could not extract notebook query constants: {', '.join(missing)}")

    chemotherapy_query = strings["CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY"]
    try:
        chemotherapy_therapy_query = chemotherapy_query.split("\nAND\n", 1)[1]
    except IndexError as error:
        raise RuntimeError("Could not separate the chemotherapy therapy query") from error

    treatment_queries = {
        "Chemotherapy": chemotherapy_therapy_query,
        "Targeted/small-molecule therapy": strings["SMALL_MOLECULE_TARGETED_THERAPY_QUERY"],
        "Antibody/ADC therapy": strings["BIOLOGIC_ANTIBODY_THERAPY_QUERY"],
        "Immunotherapy": strings["IMMUNOTHERAPY_QUERY"],
    }
    return (
        strings["PHASE_III_CNS_METASTASIS_CONTEXT"],
        treatment_queries,
        strings["PHASE_III_CRITERIA"],
        strings["PHASE_III_ARTICLE_TYPE_CRITERIA"],
        strings["TREATMENT_PUBLICATION_EXCLUSION"],
    )


def extract_phase_iii_analysis_contexts() -> dict[str, str]:
    """Return the three disease contexts used in the Phase III figures."""
    with NOTEBOOK_PATH.open(encoding="utf-8") as handle:
        notebook = json.load(handle)

    strings: dict[str, str] = {}
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        try:
            tree = ast.parse("".join(cell.get("source", [])))
        except SyntaxError:
            continue
        for statement in tree.body:
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
                continue
            target = statement.targets[0]
            if not isinstance(target, ast.Name):
                continue
            value = _static_string(statement.value)
            if value is not None:
                strings[target.id] = value

    required = {
        "CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY",
        "PRIMARY_CNS_DISEASE_CONTEXT",
        "PHASE_III_CNS_METASTASIS_CONTEXT",
    }
    missing = sorted(required.difference(strings))
    if missing:
        raise RuntimeError(f"Could not extract Phase III contexts: {', '.join(missing)}")

    overall_context = strings["CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY"].split("\nAND\n", 1)[0]
    return {
        "Aggregated neuro-oncology": overall_context,
        "Primary CNS malignancy": strings["PRIMARY_CNS_DISEASE_CONTEXT"],
        "CNS metastasis": strings["PHASE_III_CNS_METASTASIS_CONTEXT"],
    }


def extract_general_analysis_contexts() -> dict[str, str]:
    """Return the disease contexts used in the three non-Phase-III modality figures."""
    with NOTEBOOK_PATH.open(encoding="utf-8") as handle:
        notebook = json.load(handle)

    strings: dict[str, str] = {}
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        try:
            tree = ast.parse("".join(cell.get("source", [])))
        except SyntaxError:
            continue
        for statement in tree.body:
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1:
                continue
            target = statement.targets[0]
            if not isinstance(target, ast.Name):
                continue
            value = _static_string(statement.value)
            if value is not None:
                strings[target.id] = value

    required = {
        "CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY",
        "PRIMARY_CNS_DISEASE_CONTEXT",
        "CNS_METASTASIS_DISEASE_CONTEXT",
    }
    missing = sorted(required.difference(strings))
    if missing:
        raise RuntimeError(f"Could not extract treatment contexts: {', '.join(missing)}")

    overall_context = strings["CHEMOTHERAPY_NEURO_ONCOLOGY_QUERY"].split("\nAND\n", 1)[0]
    return {
        "Aggregated neuro-oncology": overall_context,
        "Primary CNS malignancy": strings["PRIMARY_CNS_DISEASE_CONTEXT"],
        "CNS metastasis": strings["CNS_METASTASIS_DISEASE_CONTEXT"],
    }


def normalize_pubmed_query(query: str) -> str:
    collapsed = " ".join(query.split())
    normalized: list[str] = []
    inside_quotes = False
    for index, character in enumerate(collapsed):
        if character == '"':
            inside_quotes = not inside_quotes
            normalized.append(character)
            continue
        if character == " " and not inside_quotes:
            previous_character = collapsed[index - 1] if index else ""
            next_character = collapsed[index + 1] if index + 1 < len(collapsed) else ""
            if previous_character == "(" or next_character == ")":
                continue
        normalized.append(character)
    return "".join(normalized)


def cached_annual_counts(query: str) -> dict[int, int]:
    query_key = hashlib.sha256(normalize_pubmed_query(query).encode("utf-8")).hexdigest()
    with sqlite3.connect(CACHE_PATH) as connection:
        rows = connection.execute(
            """
            SELECT year, count
            FROM annual_counts
            WHERE query_key = ? AND year BETWEEN ? AND ?
            ORDER BY year
            """,
            (query_key, START_YEAR, END_YEAR),
        ).fetchall()
    return {int(year): int(count) for year, count in rows}


def store_annual_count(query: str, year: int, count: int) -> None:
    """Store the retrieved count so the notebook figure uses this exact export result."""
    normalized_query = normalize_pubmed_query(query)
    query_key = hashlib.sha256(normalized_query.encode("utf-8")).hexdigest()
    with sqlite3.connect(CACHE_PATH) as connection:
        connection.execute(
            """
            INSERT INTO annual_counts
                (query_key, normalized_query, year, count, fetched_at)
            VALUES (?, ?, ?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT(query_key, year) DO UPDATE SET
                normalized_query = excluded.normalized_query,
                count = excluded.count,
                fetched_at = CURRENT_TIMESTAMP
            """,
            (query_key, normalized_query, year, count),
        )


def read_entrez(request_factory, label: str, attempts: int = 5):
    for attempt in range(1, attempts + 1):
        try:
            with request_factory() as handle:
                result = Entrez.read(handle)
            time.sleep(0.37)
            return result
        except (HTTPError, URLError, OSError, RuntimeError) as error:
            if attempt == attempts:
                raise RuntimeError(f"PubMed request failed: {label}") from error
            wait_seconds = min(30, 2**attempt)
            print(f"  Retrying {label} in {wait_seconds} seconds...", flush=True)
            time.sleep(wait_seconds)
    raise AssertionError("unreachable")


def retrieve_category_matches(
    context_query: str,
    treatment_queries: dict[str, str],
    phase_iii_criteria: str,
    article_type_criteria: str,
    publication_exclusion: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, str]]:
    matches: list[dict[str, object]] = []
    validation: list[dict[str, object]] = []
    complete_queries: dict[str, str] = {}

    for category, therapy_query in treatment_queries.items():
        complete_query = (
            f"{context_query}\nAND\n{therapy_query}\nAND\n{phase_iii_criteria}"
            f"\nAND\n{article_type_criteria}"
            f"\nNOT\n{publication_exclusion}"
        )
        complete_queries[category] = complete_query
        expected_counts = cached_annual_counts(complete_query)
        print(f"Retrieving {category} records...", flush=True)

        for year in range(START_YEAR, END_YEAR + 1):
            annual_query = f"({complete_query}) AND ({year}[pdat])"
            record = read_entrez(
                lambda annual_query=annual_query: Entrez.esearch(
                    db="pubmed",
                    term=annual_query,
                    retmax=10000,
                    sort="pub date",
                ),
                f"{category}, {year}",
            )
            pmids = [str(pmid) for pmid in record["IdList"]]
            retrieved_count = int(record["Count"])
            if retrieved_count != len(pmids):
                raise RuntimeError(
                    f"PubMed returned {len(pmids)} IDs for {retrieved_count} {category} records in {year}"
                )
            expected_count = expected_counts.get(year)
            store_annual_count(complete_query, year, retrieved_count)
            if expected_count is None:
                expected_count = retrieved_count
            matches.extend(
                {"Category": category, "Year": year, "PMID": pmid}
                for pmid in pmids
            )
            validation.append(
                {
                    "Category": category,
                    "Year": year,
                    "Figure count": expected_count if expected_count is not None else "",
                    "Retrieved records": retrieved_count,
                    "Matches figure": (
                        "Yes" if expected_count == retrieved_count else "No"
                    ) if expected_count is not None else "Not cached",
                }
            )
        print(f"  Retrieved {sum(1 for row in matches if row['Category'] == category)} matches.", flush=True)

    return matches, validation, complete_queries


def refresh_context_counts(
    context_name: str,
    context_query: str,
    treatment_queries: dict[str, str],
    required_criteria: list[str],
    excluded_criteria: list[str],
) -> None:
    """Refresh annual count-only cache entries for a treatment-modality context."""
    print(f"Refreshing {context_name} treatment-modality counts...", flush=True)
    for category, therapy_query in treatment_queries.items():
        complete_query = "\nAND\n".join([context_query, therapy_query, *required_criteria])
        for excluded_criterion in excluded_criteria:
            complete_query += f"\nNOT\n{excluded_criterion}"
        category_total = 0
        print(f"  Retrieving {category} counts...", flush=True)
        for year in range(START_YEAR, END_YEAR + 1):
            annual_query = f"({complete_query}) AND ({year}[pdat])"
            record = read_entrez(
                lambda annual_query=annual_query: Entrez.esearch(
                    db="pubmed",
                    term=annual_query,
                    retmax=0,
                ),
                f"{context_name}, {category}, {year}",
            )
            count = int(record["Count"])
            store_annual_count(complete_query, year, count)
            category_total += count
        print(f"    Cached {category_total} publications.", flush=True)


def clean_title(title: object) -> str:
    value = html.unescape(str(title))
    value = re.sub(r"<[^>]+>", "", value)
    return " ".join(value.split())


def fetch_titles(pmids: list[str]) -> dict[str, str]:
    titles: dict[str, str] = {}
    batch_size = 200
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start : start + batch_size]
        print(
            f"Fetching titles {start + 1}-{min(start + batch_size, len(pmids))} of {len(pmids)}...",
            flush=True,
        )
        summaries = read_entrez(
            lambda batch=batch: Entrez.esummary(
                db="pubmed", id=",".join(batch), retmode="xml"
            ),
            f"titles for PMID batch starting at {start + 1}",
        )
        for summary in summaries:
            pmid = str(summary.attributes.get("uid") or summary.get("Id", ""))
            if pmid:
                titles[pmid] = clean_title(summary.get("Title", ""))
    missing = sorted(set(pmids).difference(titles))
    if missing:
        raise RuntimeError(f"Missing titles for {len(missing)} PMID(s): {', '.join(missing[:10])}")
    return titles


def column_name(index: int) -> str:
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result


def clean_xml_text(value: object) -> str:
    text = str(value)
    text = re.sub(r"[\x00-\x08\x0B\x0C\x0E-\x1F]", "", text)
    return text[:32767]


def cell_xml(row: int, column: int, value: object, style: int = 0) -> str:
    reference = f"{column_name(column)}{row}"
    style_attribute = f' s="{style}"' if style else ""
    if isinstance(value, int):
        return f'<c r="{reference}"{style_attribute}><v>{value}</v></c>'
    text = clean_xml_text(value)
    return (
        f'<c r="{reference}"{style_attribute} t="inlineStr">'
        f'<is><t xml:space="preserve">{escape(text)}</t></is></c>'
    )


def worksheet_xml(
    rows: list[list[object]],
    widths: list[float],
    hyperlink_column: int | None = None,
    wrap_columns: set[int] | None = None,
) -> tuple[str, str | None]:
    wrap_columns = wrap_columns or set()
    row_xml: list[str] = []
    hyperlink_xml: list[str] = []
    relationship_xml: list[str] = []
    relationship_index = 1

    for row_index, row in enumerate(rows, start=1):
        cells: list[str] = []
        for column_index, value in enumerate(row, start=1):
            if row_index == 1:
                style = 1
            elif hyperlink_column == column_index:
                style = 2
            elif column_index in wrap_columns:
                style = 3
            else:
                style = 0
            cells.append(cell_xml(row_index, column_index, value, style))
        row_xml.append(f'<row r="{row_index}">{"".join(cells)}</row>')

        if hyperlink_column and row_index > 1 and len(row) >= hyperlink_column:
            target = clean_xml_text(row[hyperlink_column - 1])
            if target:
                reference = f"{column_name(hyperlink_column)}{row_index}"
                relationship_id = f"rId{relationship_index}"
                hyperlink_xml.append(f'<hyperlink ref="{reference}" r:id="{relationship_id}"/>')
                relationship_xml.append(
                    '<Relationship '
                    f'Id="{relationship_id}" '
                    'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" '
                    f'Target={quoteattr(target)} TargetMode="External"/>'
                )
                relationship_index += 1

    last_column = column_name(max(len(row) for row in rows))
    last_row = len(rows)
    columns = "".join(
        f'<col min="{index}" max="{index}" width="{width}" customWidth="1"/>'
        for index, width in enumerate(widths, start=1)
    )
    hyperlinks = f'<hyperlinks>{"".join(hyperlink_xml)}</hyperlinks>' if hyperlink_xml else ""
    worksheet = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="A1:{last_column}{last_row}"/>'
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        '</sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="15"/>'
        f'<cols>{columns}</cols>'
        f'<sheetData>{"".join(row_xml)}</sheetData>'
        f'<autoFilter ref="A1:{last_column}{last_row}"/>'
        f'{hyperlinks}'
        '</worksheet>'
    )
    relationships = None
    if relationship_xml:
        relationships = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{"".join(relationship_xml)}'
            '</Relationships>'
        )
    return worksheet, relationships


def write_workbook(
    output_path: Path,
    studies: list[dict[str, object]],
    category_matches: list[dict[str, object]],
    validation: list[dict[str, object]],
    complete_queries: dict[str, str],
) -> None:
    studies_rows: list[list[object]] = [["Year", "Title", "URL"]]
    studies_rows.extend([[row["Year"], row["Title"], row["URL"]] for row in studies])

    match_rows: list[list[object]] = [["Treatment category", "Year", "PMID", "Title", "URL"]]
    match_rows.extend(
        [[row["Category"], row["Year"], row["PMID"], row["Title"], row["URL"]]
         for row in category_matches]
    )

    validation_rows: list[list[object]] = [
        ["Treatment category", "Year", "Figure count", "Retrieved records", "Matches figure"]
    ]
    validation_rows.extend(
        [[row["Category"], row["Year"], row["Figure count"], row["Retrieved records"], row["Matches figure"]]
         for row in validation]
    )

    query_rows: list[list[object]] = [["Treatment category", "Complete PubMed query"]]
    query_rows.extend([[category, query] for category, query in complete_queries.items()])

    sheets = [
        ("Studies", studies_rows, [10, 95, 42], 3, {2}),
        ("Category matches", match_rows, [34, 10, 14, 95, 42], 5, {4}),
        ("Count validation", validation_rows, [34, 10, 15, 18, 16], None, set()),
        ("PubMed queries", query_rows, [34, 120], None, {2}),
    ]

    content_type_overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        f'{content_type_overrides}</Types>'
    )
    root_relationships = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
        'Target="docProps/app.xml"/>'
        '</Relationships>'
    )
    sheet_entries = "".join(
        f'<sheet name={quoteattr(name)} sheetId="{index}" r:id="rId{index}"/>'
        for index, (name, *_rest) in enumerate(sheets, start=1)
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<bookViews><workbookView/></bookViews>'
        f'<sheets>{sheet_entries}</sheets>'
        '<calcPr calcId="191029" fullCalcOnLoad="1"/></workbook>'
    )
    workbook_relationships = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(
            '<Relationship '
            f'Id="rId{index}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{index}.xml"/>'
            for index in range(1, len(sheets) + 1)
        )
        + f'<Relationship Id="rId{len(sheets) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        '</Relationships>'
    )
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="3">'
        '<font><sz val="11"/><name val="Calibri"/><family val="2"/></font>'
        '<font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font>'
        '<font><u/><color rgb="FF0563C1"/><sz val="11"/><name val="Calibri"/></font>'
        '</fonts>'
        '<fills count="3"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/>'
        '<bgColor indexed="64"/></patternFill></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="4">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
        '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
        '<alignment vertical="top" wrapText="1"/></xf>'
        '</cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )
    timestamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    core_properties = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties '
        'xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        '<dc:title>Phase III CNS metastasis studies, 2000-2025</dc:title>'
        '<dc:creator>PubMed export</dc:creator>'
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{timestamp}</dcterms:created>'
        f'<dcterms:modified xsi:type="dcterms:W3CDTF">{timestamp}</dcterms:modified>'
        '</cp:coreProperties>'
    )
    app_properties = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        '<Application>Microsoft Excel Compatible</Application>'
        f'<TitlesOfParts><vt:vector size="{len(sheets)}" baseType="lpstr">'
        + "".join(f'<vt:lpstr>{escape(name)}</vt:lpstr>' for name, *_rest in sheets)
        + '</vt:vector></TitlesOfParts></Properties>'
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as workbook_zip:
        workbook_zip.writestr("[Content_Types].xml", content_types)
        workbook_zip.writestr("_rels/.rels", root_relationships)
        workbook_zip.writestr("docProps/core.xml", core_properties)
        workbook_zip.writestr("docProps/app.xml", app_properties)
        workbook_zip.writestr("xl/workbook.xml", workbook)
        workbook_zip.writestr("xl/_rels/workbook.xml.rels", workbook_relationships)
        workbook_zip.writestr("xl/styles.xml", styles)
        for index, (_name, rows, widths, hyperlink_column, wrap_columns) in enumerate(sheets, start=1):
            sheet, relationships = worksheet_xml(rows, widths, hyperlink_column, wrap_columns)
            workbook_zip.writestr(f"xl/worksheets/sheet{index}.xml", sheet)
            if relationships:
                workbook_zip.writestr(
                    f"xl/worksheets/_rels/sheet{index}.xml.rels",
                    relationships,
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", "user@example.com"),
        help="Email supplied to NCBI Entrez (default: NCBI_EMAIL or notebook placeholder)",
    )
    parser.add_argument(
        "--refresh-all-contexts",
        action="store_true",
        help="Refresh cached counts for aggregated, primary-CNS, and CNS-metastasis Phase III figures",
    )
    parser.add_argument(
        "--refresh-general-contexts",
        action="store_true",
        help="Refresh cached counts for the three non-Phase-III treatment-modality figures",
    )
    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.tool = "lit_trends_phase_iii_cns_metastasis_export"
    Entrez.max_tries = 1

    (
        context_query,
        treatment_queries,
        phase_iii_criteria,
        article_type_criteria,
        publication_exclusion,
    ) = extract_notebook_queries()
    if args.refresh_general_contexts:
        general_contexts = extract_general_analysis_contexts()
        for context_name, general_context in general_contexts.items():
            refresh_context_counts(
                context_name,
                general_context,
                treatment_queries,
                [],
                [publication_exclusion],
            )
    if args.refresh_all_contexts:
        contexts = extract_phase_iii_analysis_contexts()
        for context_name in ("Aggregated neuro-oncology", "Primary CNS malignancy"):
            print(f"Refreshing {context_name} Phase III counts...", flush=True)
            retrieve_category_matches(
                contexts[context_name],
                treatment_queries,
                phase_iii_criteria,
                article_type_criteria,
                publication_exclusion,
            )
    matches, validation, complete_queries = retrieve_category_matches(
        context_query,
        treatment_queries,
        phase_iii_criteria,
        article_type_criteria,
        publication_exclusion,
    )

    unique_pmids = sorted({str(row["PMID"]) for row in matches}, key=int)
    titles = fetch_titles(unique_pmids)
    for row in matches:
        pmid = str(row["PMID"])
        row["Title"] = titles[pmid]
        row["URL"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

    category_order = {category: index for index, category in enumerate(treatment_queries)}
    matches.sort(
        key=lambda row: (
            int(row["Year"]),
            category_order[str(row["Category"])],
            str(row["Title"]).casefold(),
            int(str(row["PMID"])),
        )
    )

    study_by_pmid: dict[str, dict[str, object]] = {}
    for row in matches:
        pmid = str(row["PMID"])
        study_by_pmid.setdefault(
            pmid,
            {
                "Year": row["Year"],
                "Title": row["Title"],
                "URL": row["URL"],
            },
        )
    studies = sorted(
        study_by_pmid.values(),
        key=lambda row: (int(row["Year"]), str(row["Title"]).casefold()),
    )

    write_workbook(args.output, studies, matches, validation, complete_queries)
    mismatches = [row for row in validation if row["Matches figure"] == "No"]
    print(f"Created {args.output}", flush=True)
    print(f"Unique studies: {len(studies)}", flush=True)
    print(f"Category-study matches: {len(matches)}", flush=True)
    print(f"Category-year count mismatches versus figure cache: {len(mismatches)}", flush=True)


if __name__ == "__main__":
    main()
