# TCGA-BRCA Hot/Cold immune-TME preparation

## Data provenance

- Source: NCI Genomic Data Commons (GDC), official API
- Files endpoint: https://api.gdc.cancer.gov/files
- Data endpoint: https://api.gdc.cancer.gov/data/{file_uuid}
- Project: TCGA-BRCA
- Data category/type: Transcriptome Profiling / Gene Expression Quantification
- Experimental strategy: RNA-Seq
- Workflow: STAR - Counts
- Sample type: Primary Tumor
- Access: open
- GDC release/status: Data Release 46.0 - August 10, 2026 (recorded during GDC metadata validation)
- Retrieval date (UTC): 2026-09-30T18:43:58+00:00

## Retrieval and processing

The exact API request and response are saved in `data/processed/gdc_brca_star_counts_api_query.json` and `data/processed/gdc_brca_star_counts_api_response.json`. The file-level manifest is `data/processed/gdc_brca_star_counts_manifest.csv`. Raw files are named by GDC file UUID and were downloaded only from the GDC `/data/{file_uuid}` endpoint.

The metadata query returned 1111 qualifying files; 1106 files were selected after duplicate-sample handling and 1106 were downloaded during this run. The STAR-Counts schema was inspected from a downloaded file before parsing. Its resolved columns were: gene identifier `gene_id`, gene symbol `gene_name`, and unstranded counts `unstranded`.

The processed matrix has 27 found panel genes by 1106 unique TCGA samples. Missing genes: None.

`BRCA_27gene_raw_counts.csv` contains raw unstranded gene-level counts. `BRCA_27gene_log2_counts.csv` contains `log2(counts + 1)` transformed counts only. It is not TPM, FPKM, or FPKM-UQ.

Summary rows are removed by retaining Ensembl gene IDs beginning `ENSG`. Where multiple retained Ensembl IDs share a gene symbol, their unstranded counts are summed per sample; this is an explicit aggregation, not a silent row drop. Duplicate samples are recorded in the manifest. The deterministic rule retains the newest `updated_datetime`, with UUID as a tie-breaker.

## Outputs

- `scripts/download_tcga_paad.py`: reproducible query, validation, download, and matrix build script
- `data/processed/BRCA_27gene_raw_counts.csv`: raw count matrix
- `data/processed/BRCA_27gene_log2_counts.csv`: log2-transformed count matrix
- `data/processed/BRCA_sample_metadata.csv`: sample and GDC file metadata
