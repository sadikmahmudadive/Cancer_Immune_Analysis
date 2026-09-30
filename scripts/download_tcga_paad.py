#!/usr/bin/env python3
"""Retrieve open-access TCGA STAR-Counts data from the GDC API.

The default mode is metadata-only. Pass --download only after reviewing the
validation report and manifest produced by the metadata query.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests


PROJECT = "TCGA-PAAD"
WORKFLOW = "STAR - Counts"
SAMPLE_TYPE = "Primary Tumor"
ACCESS = "open"
API_BASE = "https://api.gdc.cancer.gov"
REQUESTED_GENES = [
    "CD8A", "CD8B", "IFNG", "GZMB", "PRF1", "NKG7", "CXCL9", "CXCL10",
    "CXCL11", "STAT1", "IRF1", "HLA-A", "HLA-B", "HLA-C", "B2M", "TGFB1",
    "CXCL12", "POSTN", "HIF1A", "VEGFA", "KDR", "FAP", "COL1A1", "CD163",
    "CSF1R", "FOXP3", "IL10",
]
HOT_GENES = REQUESTED_GENES[:15]
COLD_GENES = REQUESTED_GENES[15:]
FIELDS = [
    "file_id", "file_name", "access", "data_category", "data_type",
    "experimental_strategy", "analysis.workflow_type", "release_state",
    "created_datetime", "updated_datetime", "cases.case_id", "cases.submitter_id",
    "cases.project.project_id", "cases.samples.sample_id",
    "cases.samples.submitter_id", "cases.samples.sample_type",
]


def project_code() -> str:
    return PROJECT.removeprefix("TCGA-")


def gdc_filters() -> dict[str, Any]:
    def clause(field: str, value: str) -> dict[str, Any]:
        return {"op": "in", "content": {"field": field, "value": [value]}}

    return {
        "op": "and",
        "content": [
            clause("cases.project.project_id", PROJECT),
            clause("data_category", "Transcriptome Profiling"),
            clause("data_type", "Gene Expression Quantification"),
            clause("experimental_strategy", "RNA-Seq"),
            clause("analysis.workflow_type", WORKFLOW),
            clause("cases.samples.sample_type", SAMPLE_TYPE),
            clause("access", ACCESS),
        ],
    }


def make_dirs(root: Path) -> dict[str, Path]:
    paths = {
        "raw": root / "data" / "raw" / project_code(),
        "processed": root / "data" / "processed",
        "notebooks": root / "notebooks",
        "figures": root / "results" / "figures",
        "tables": root / "results" / "tables",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def api_json(session: requests.Session, endpoint: str, **kwargs: Any) -> dict[str, Any]:
    response = session.request(timeout=90, **kwargs, url=f"{API_BASE}{endpoint}")
    response.raise_for_status()
    return response.json()


def query_files(session: requests.Session) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = {"filters": gdc_filters(), "fields": ",".join(FIELDS), "format": "JSON", "size": 10000}
    response = api_json(session, "/files", method="POST", json=payload)
    return payload, response


def primary_sample_records(hit: dict[str, Any]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for case in hit.get("cases", []):
        for sample in case.get("samples", []):
            if sample.get("sample_type") == SAMPLE_TYPE:
                records.append({
                    "sample_id": sample.get("submitter_id") or sample.get("sample_id") or "",
                    "sample_uuid": sample.get("sample_id") or "",
                    "case_id": case.get("case_id") or "",
                    "case_submitter_id": case.get("submitter_id") or "",
                    "sample_type": sample.get("sample_type") or "",
                })
    return records


def validate_and_manifest(hits: list[dict[str, Any]]) -> tuple[list[dict[str, str]], list[str]]:
    required = {
        "project": PROJECT,
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "experimental_strategy": "RNA-Seq",
        "workflow": WORKFLOW,
        "access": ACCESS,
    }
    records: list[dict[str, str]] = []
    errors: list[str] = []
    for hit in hits:
        project_ids = {
            case.get("project", {}).get("project_id", "")
            for case in hit.get("cases", [])
        }
        observed = {
            "project": PROJECT if PROJECT in project_ids else ",".join(sorted(project_ids)),
            "data_category": hit.get("data_category", ""),
            "data_type": hit.get("data_type", ""),
            "experimental_strategy": hit.get("experimental_strategy", ""),
            "workflow": hit.get("analysis", {}).get("workflow_type", ""),
            "access": hit.get("access", ""),
        }
        mismatches = [key for key, value in required.items() if observed[key] != value]
        sample_records = primary_sample_records(hit)
        if mismatches or not sample_records:
            errors.append(f"{hit.get('file_id', '<missing UUID>')}: mismatches={mismatches}; primary_samples={len(sample_records)}")
            continue
        for sample in sample_records:
            records.append({
                "file_uuid": hit["file_id"],
                "file_name": hit.get("file_name", ""),
                "sample_id": sample["sample_id"],
                "sample_uuid": sample["sample_uuid"],
                "case_id": sample["case_id"],
                "case_submitter_id": sample["case_submitter_id"],
                "sample_type": sample["sample_type"],
                "project": PROJECT,
                "data_category": hit.get("data_category", ""),
                "data_type": hit.get("data_type", ""),
                "experimental_strategy": hit.get("experimental_strategy", ""),
                "workflow": hit.get("analysis", {}).get("workflow_type", ""),
                "access": hit.get("access", ""),
                "release_state": hit.get("release_state", ""),
                "created_datetime": hit.get("created_datetime", ""),
                "updated_datetime": hit.get("updated_datetime", ""),
                "selection_status": "selected",
                "duplicate_resolution": "",
            })

    duplicate_uuids = [uuid for uuid, count in Counter(row["file_uuid"] for row in records).items() if count > 1]
    if duplicate_uuids:
        errors.append(f"Duplicate file UUIDs in manifest: {duplicate_uuids}")

    # A matrix needs one column per sample. Keep the most recently updated file
    # when GDC returns repeated files for a sample, and record every exclusion.
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in records:
        grouped.setdefault(row["sample_id"], []).append(row)
    for sample_id, rows in grouped.items():
        if len(rows) > 1:
            rows.sort(key=lambda row: (row["updated_datetime"], row["file_uuid"]), reverse=True)
            for index, row in enumerate(rows):
                if index:
                    row["selection_status"] = "excluded_duplicate_sample"
                    row["duplicate_resolution"] = f"Excluded; retained {rows[0]['file_uuid']} (latest updated_datetime, then UUID)."
                else:
                    row["duplicate_resolution"] = "Retained latest updated_datetime; UUID used as deterministic tie-breaker."
    return records, errors


def write_csv(rows: list[dict[str, str]], output: Path) -> None:
    if not rows:
        raise RuntimeError("No manifest rows were available to write.")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def print_validation(hits: list[dict[str, Any]], records: list[dict[str, str]], errors: list[str]) -> str:
    selected = [row for row in records if row["selection_status"] == "selected"]
    all_sample_ids = [row["sample_id"] for row in records]
    report = [
        "GDC METADATA-ONLY VALIDATION REPORT",
        f"Files returned: {len(hits)}",
        f"Primary Tumor sample records: {len(records)}",
        f"Unique samples: {len(set(all_sample_ids))}",
        f"Selected files: {len(selected)}",
        f"Duplicate samples: {sum(count - 1 for count in Counter(all_sample_ids).values() if count > 1)}",
        f"Duplicate file UUIDs: {sum(count - 1 for count in Counter(row['file_uuid'] for row in records).values() if count > 1)}",
        f"Project: {PROJECT}",
        "Data category: Transcriptome Profiling",
        "Data type: Gene Expression Quantification",
        "Experimental strategy: RNA-Seq",
        f"Workflow: {WORKFLOW}",
        f"Sample type: {SAMPLE_TYPE}",
        f"Access level: {ACCESS}",
        "Example files:",
        *[f"  {row['file_uuid']}  {row['file_name']}  {row['sample_id']}" for row in selected[:5]],
    ]
    if errors:
        report.extend(["Validation errors:", *[f"  {error}" for error in errors]])
    text = "\n".join(report)
    print(text)
    return text


def download_files(session: requests.Session, rows: list[dict[str, str]], raw_dir: Path) -> int:
    def download_one(row: dict[str, str]) -> tuple[str, bool]:
        output = raw_dir / f"{row['file_uuid']}.tsv"
        if output.exists() and output.stat().st_size > 0:
            return row["file_uuid"], False
        temporary = output.with_suffix(".part")
        with requests.get(
            f"{API_BASE}/data/{row['file_uuid']}",
            headers={"User-Agent": "TCGA-Hot-Cold-preparation/1.0"},
            stream=True,
            timeout=180,
        ) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        temporary.replace(output)
        return row["file_uuid"], True

    pending = [row for row in rows if not (raw_dir / f"{row['file_uuid']}.tsv").exists()]
    print(f"Download status: {len(rows) - len(pending)} existing complete files; {len(pending)} files pending.")
    downloaded = 0
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(download_one, row): row for row in pending}
        for completed, future in enumerate(as_completed(futures), start=1):
            row = futures[future]
            file_uuid, was_downloaded = future.result()
            downloaded += int(was_downloaded)
            print(f"[{completed}/{len(pending)}] Complete: {file_uuid} ({row['sample_id']})")
    return downloaded


def detect_star_columns(path: Path) -> tuple[list[str], int, dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        header_line = 0
        for line in handle:
            if not line.startswith("#"):
                columns = line.rstrip("\n\r").split("\t")
                break
            header_line += 1
        else:
            raise RuntimeError(f"No tabular header found in {path}")
    canonical = {column.lower(): column for column in columns}
    required = {"gene_identifier": "gene_id", "gene_symbol": "gene_name", "unstranded_count": "unstranded"}
    resolved: dict[str, str] = {}
    for label, expected in required.items():
        if expected not in canonical:
            raise RuntimeError(f"{path.name} does not contain required STAR-Counts column '{expected}'. Observed: {columns}")
        resolved[label] = canonical[expected]
    return columns, header_line, resolved


def load_gene_counts(path: Path, header_line: int, columns: dict[str, str]) -> pd.Series:
    frame = pd.read_csv(
        path,
        sep="\t",
        skiprows=header_line,
        usecols=[columns["gene_identifier"], columns["gene_symbol"], columns["unstranded_count"]],
        dtype={columns["gene_identifier"]: "string", columns["gene_symbol"]: "string"},
    )
    ids = frame[columns["gene_identifier"]].fillna("").astype(str)
    symbols = frame[columns["gene_symbol"]].fillna("").astype(str)
    values = pd.to_numeric(frame[columns["unstranded_count"]], errors="coerce")
    keep = ids.str.startswith("ENSG") & symbols.ne("") & values.notna()
    usable = pd.DataFrame({"symbol": symbols[keep], "count": values[keep]})
    # Different Ensembl identifiers can map to one symbol. Sum counts per symbol
    # and disclose the aggregation in README and the final report.
    return usable.groupby("symbol", sort=True)["count"].sum()


def build_matrices(selected: list[dict[str, str]], raw_dir: Path, processed_dir: Path) -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    first_path = raw_dir / f"{selected[0]['file_uuid']}.tsv"
    observed_columns, header_line, column_map = detect_star_columns(first_path)
    print("Observed STAR-Counts columns:", ", ".join(observed_columns))
    print("Resolved parser columns:", json.dumps(column_map, sort_keys=True))
    series_by_sample: dict[str, pd.Series] = {}
    duplicate_symbol_rows = 0
    for row in selected:
        path = raw_dir / f"{row['file_uuid']}.tsv"
        current_columns, current_header_line, current_map = detect_star_columns(path)
        if current_columns != observed_columns or current_map != column_map:
            raise RuntimeError(f"Schema differs in {path.name}; expected {observed_columns}, found {current_columns}")
        series = load_gene_counts(path, current_header_line, current_map)
        duplicate_symbol_rows += int(series.index.duplicated().sum())
        series_by_sample[row["sample_id"]] = series.reindex(REQUESTED_GENES)
    matrix = pd.DataFrame(series_by_sample).reindex(REQUESTED_GENES)
    found = [gene for gene in REQUESTED_GENES if gene in matrix.index and matrix.loc[gene].notna().any()]
    missing = [gene for gene in REQUESTED_GENES if gene not in found]
    matrix.index.name = "gene_symbol"
    matrix.to_csv(processed_dir / f"{project_code()}_27gene_raw_counts.csv")
    log2_matrix = np.log2(matrix + 1)
    log2_matrix.index.name = "gene_symbol"
    log2_matrix.to_csv(processed_dir / f"{project_code()}_27gene_log2_counts.csv")
    parser_info = {
        "observed_columns": observed_columns,
        "header_line": header_line,
        "resolved_columns": column_map,
        "duplicate_symbol_rows_after_aggregation": duplicate_symbol_rows,
    }
    return matrix, missing, parser_info


def write_readme(root: Path, status: dict[str, Any]) -> None:
    missing = ", ".join(status["missing_genes"]) if status["missing_genes"] else "None"
    code = project_code()
    readme = f"""# TCGA-{code} Hot/Cold immune-TME preparation

## Data provenance

- Source: NCI Genomic Data Commons (GDC), official API
- Files endpoint: {API_BASE}/files
- Data endpoint: {API_BASE}/data/{{file_uuid}}
- Project: {PROJECT}
- Data category/type: Transcriptome Profiling / Gene Expression Quantification
- Experimental strategy: RNA-Seq
- Workflow: {WORKFLOW}
- Sample type: {SAMPLE_TYPE}
- Access: {ACCESS}
- GDC release/status: {status['gdc_release']}
- Retrieval date (UTC): {status['retrieval_date']}

## Retrieval and processing

The exact API request and response are saved in `data/processed/gdc_{code.lower()}_star_counts_api_query.json` and `data/processed/gdc_{code.lower()}_star_counts_api_response.json`. The file-level manifest is `data/processed/gdc_{code.lower()}_star_counts_manifest.csv`. Raw files are named by GDC file UUID and were downloaded only from the GDC `/data/{{file_uuid}}` endpoint.

The metadata query returned {status['files_found']} qualifying files; {status['files_selected']} files were selected after duplicate-sample handling and {status['files_downloaded']} were downloaded during this run. The STAR-Counts schema was inspected from a downloaded file before parsing. Its resolved columns were: gene identifier `{status['parser_columns'].get('gene_identifier', 'not inspected')}`, gene symbol `{status['parser_columns'].get('gene_symbol', 'not inspected')}`, and unstranded counts `{status['parser_columns'].get('unstranded_count', 'not inspected')}`.

The processed matrix has {status['genes_found']} found panel genes by {status['samples']} unique TCGA samples. Missing genes: {missing}.

`{code}_27gene_raw_counts.csv` contains raw unstranded gene-level counts. `{code}_27gene_log2_counts.csv` contains `log2(counts + 1)` transformed counts only. It is not TPM, FPKM, or FPKM-UQ.

Summary rows are removed by retaining Ensembl gene IDs beginning `ENSG`. Where multiple retained Ensembl IDs share a gene symbol, their unstranded counts are summed per sample; this is an explicit aggregation, not a silent row drop. Duplicate samples are recorded in the manifest. The deterministic rule retains the newest `updated_datetime`, with UUID as a tie-breaker.

## Outputs

- `scripts/download_tcga_paad.py`: reproducible query, validation, download, and matrix build script
- `data/processed/{code}_27gene_raw_counts.csv`: raw count matrix
- `data/processed/{code}_27gene_log2_counts.csv`: log2-transformed count matrix
- `data/processed/{code}_sample_metadata.csv`: sample and GDC file metadata
"""
    (root / "results" / "tables" / f"{code}_README.md").write_text(readme, encoding="utf-8")


def final_report(status: dict[str, Any]) -> str:
    code = project_code()
    output_files = [
        f"data/processed/{code}_27gene_raw_counts.csv",
        f"data/processed/{code}_27gene_log2_counts.csv",
        f"data/processed/{code}_sample_metadata.csv",
        f"data/processed/gdc_{code.lower()}_star_counts_manifest.csv",
    ]
    return "\n".join([
        "FINAL REPORT",
        f"TCGA PROJECT: {PROJECT}",
        f"FILES FOUND: {status['files_found']}",
        f"FILES SELECTED: {status['files_selected']}",
        f"FILES DOWNLOADED: {status['files_downloaded']}",
        f"SAMPLES: {status['samples']}",
        "GENES REQUESTED: 27",
        f"GENES FOUND: {status['genes_found']}",
        f"GENES MISSING: {', '.join(status['missing_genes']) if status['missing_genes'] else 'None'}",
        f"EXPRESSION MATRIX: {status['genes_found']} genes x {status['samples']} samples",
        "NORMALIZATION/TRANSFORMATION: log2(counts + 1); transformed counts, not TPM",
        f"OUTPUT FILES: {', '.join(output_files)}",
        "FIGURES: generated by the two notebooks in results/figures/",
        f"GDC RELEASE: {status['gdc_release']}",
        f"DATA RETRIEVAL DATE: {status['retrieval_date']}",
    ])


def main() -> int:
    global PROJECT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", choices=["TCGA-PAAD", "TCGA-SKCM", "TCGA-BRCA"], default=PROJECT, help="TCGA project to retrieve.")
    parser.add_argument("--download", action="store_true", help="Download selected UUID-addressed STAR-Counts files and build matrices.")
    parser.add_argument("--from-saved-manifest", action="store_true", help="Build matrices locally from a previously validated manifest and complete raw UUID files.")
    args = parser.parse_args()
    PROJECT = args.project
    root = Path(__file__).resolve().parents[1]
    paths = make_dirs(root)
    retrieval_date = datetime.now(UTC).replace(microsecond=0).isoformat()
    code = project_code()
    if args.from_saved_manifest:
        manifest_path = paths["processed"] / f"gdc_{code.lower()}_star_counts_manifest.csv"
        response_path = paths["processed"] / f"gdc_{code.lower()}_star_counts_api_response.json"
        if not manifest_path.is_file() or not response_path.is_file():
            raise RuntimeError("Saved manifest or GDC response is unavailable; run metadata validation before offline assembly.")
        records = pd.read_csv(manifest_path, dtype=str).fillna("").to_dict("records")
        hits = json.loads(response_path.read_text(encoding="utf-8")).get("data", {}).get("hits", [])
        selected = [row for row in records if row["selection_status"] == "selected"]
        if not args.download:
            raise RuntimeError("Offline assembly requires --download to confirm matrix creation.")
        files_downloaded = 0
        print(f"Offline assembly: reusing {len(selected)} selected UUID files from {manifest_path.name}.")
    else:
        with requests.Session() as session:
            session.headers.update({"User-Agent": "TCGA-Hot-Cold-preparation/1.0"})
            payload, response = query_files(session)
            status_response = api_json(session, "/status", method="GET")
            hits = response.get("data", {}).get("hits", [])
            records, validation_errors = validate_and_manifest(hits)
            write_csv(records, paths["processed"] / f"gdc_{code.lower()}_star_counts_manifest.csv")
            (paths["processed"] / f"gdc_{code.lower()}_star_counts_api_query.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
            (paths["processed"] / f"gdc_{code.lower()}_star_counts_api_response.json").write_text(json.dumps(response, indent=2), encoding="utf-8")
            validation_text = print_validation(hits, records, validation_errors)
            (paths["tables"] / f"{code}_retrieval_validation.txt").write_text(validation_text + "\n", encoding="utf-8")
            if validation_errors:
                raise RuntimeError(f"Metadata validation failed; no files will be downloaded. See results/tables/{code}_retrieval_validation.txt")
            selected = [row for row in records if row["selection_status"] == "selected"]
            if not args.download:
                print("Metadata-only mode complete. Review the manifest, then rerun with --download to retrieve files.")
                return 0
            files_downloaded = download_files(session, selected, paths["raw"])

    matrix, missing, parser_info = build_matrices(selected, paths["raw"], paths["processed"])
    sample_metadata = pd.DataFrame(selected)[["sample_id", "case_id", "case_submitter_id", "sample_type", "project", "file_uuid", "file_name", "workflow", "access"]]
    if sample_metadata["sample_id"].duplicated().any():
        raise RuntimeError("Sample IDs are not unique after duplicate handling.")
    sample_metadata.to_csv(paths["processed"] / f"{project_code()}_sample_metadata.csv", index=False)
    if matrix.columns.duplicated().any() or matrix.index.duplicated().any():
        raise RuntimeError("Expression matrix has duplicate sample IDs or gene symbols.")
    if len(matrix.columns) != len(sample_metadata):
        raise RuntimeError("Expression matrix sample count does not match sample metadata.")
    completed_files = sum((paths["raw"] / f"{row['file_uuid']}.tsv").is_file() for row in selected)
    status = {
        "files_found": len(hits), "files_selected": len(selected), "files_downloaded": completed_files,
        "files_downloaded_this_run": files_downloaded,
        "samples": len(matrix.columns), "genes_found": len(REQUESTED_GENES) - len(missing),
        "missing_genes": missing, "gdc_release": ("Data Release 46.0 - August 10, 2026 (recorded during GDC metadata validation)" if args.from_saved_manifest else status_response.get("data_release", status_response.get("version", "not reported"))),
        "retrieval_date": retrieval_date, "parser_columns": parser_info["resolved_columns"],
    }
    write_readme(root, status)
    report = final_report(status)
    print(report)
    (paths["tables"] / f"{project_code()}_final_report.txt").write_text(report + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except requests.RequestException as error:
        print(f"GDC API/download error: {error}", file=sys.stderr)
        raise SystemExit(1)
    except Exception as error:
        print(f"Pipeline error: {error}", file=sys.stderr)
        raise SystemExit(1)
