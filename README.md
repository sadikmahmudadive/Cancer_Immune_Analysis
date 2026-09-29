# TCGA-PAAD Hot/Cold immune-TME preparation

## Data provenance

- Source: NCI Genomic Data Commons (GDC), official API
- Files endpoint: https://api.gdc.cancer.gov/files
- Data endpoint: https://api.gdc.cancer.gov/data/{file_uuid}
- Project: TCGA-PAAD
- Data category/type: Transcriptome Profiling / Gene Expression Quantification
- Experimental strategy: RNA-Seq
- Workflow: STAR - Counts
- Sample type: Primary Tumor
- Access: open
- GDC release/status: Data Release 46.0 - August 10, 2026
- Retrieval date (UTC): 2026-09-29T20:01:58+00:00

## Retrieval and processing

The exact API request and response are saved in `data/processed/gdc_paad_star_counts_api_query.json` and `data/processed/gdc_paad_star_counts_api_response.json`. The file-level manifest is `data/processed/gdc_paad_star_counts_manifest.csv`. Raw files are named by GDC file UUID and were downloaded only from the GDC `/data/{file_uuid}` endpoint.

The metadata query returned 178 qualifying files; 178 files were selected after duplicate-sample handling and 178 were downloaded during this run. The STAR-Counts schema was inspected from a downloaded file before parsing. Its resolved columns were: gene identifier `gene_id`, gene symbol `gene_name`, and unstranded counts `unstranded`.

The processed matrix has 27 found panel genes by 178 unique TCGA samples. Missing genes: None.

`PAAD_27gene_raw_counts.csv` contains raw unstranded gene-level counts. `PAAD_27gene_log2_counts.csv` contains `log2(counts + 1)` transformed counts only. It is not TPM, FPKM, or FPKM-UQ.

Summary rows are removed by retaining Ensembl gene IDs beginning `ENSG`. Where multiple retained Ensembl IDs share a gene symbol, their unstranded counts are summed per sample; this is an explicit aggregation, not a silent row drop. Duplicate samples are recorded in the manifest. The deterministic rule retains the newest `updated_datetime`, with UUID as a tie-breaker.

## Outputs

- `scripts/download_tcga_paad.py`: reproducible query, validation, download, and matrix build script
- `data/processed/PAAD_27gene_raw_counts.csv`: raw count matrix
- `data/processed/PAAD_27gene_log2_counts.csv`: log2-transformed count matrix
- `data/processed/PAAD_sample_metadata.csv`: sample and GDC file metadata
- `notebooks/01_PAAD_Data_Preparation.ipynb`: preparation quality checks
- `notebooks/02_PAAD_Hot_Cold_Panel.ipynb`: continuous Hot/Cold scoring without classification
