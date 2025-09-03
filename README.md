# Identification and Characterization of Novel Chikungunya Virus Polymerase Inhibitors
This repository contains the code and data used in the publication **"Identification and Characterization of Novel Chikungunya Virus Polymerase Inhibitors"** by **Yin** et al. (2025).

## What is RAVA?
RAVA is derived from [LAVA](https://www.biorxiv.org/content/10.1101/2019.12.17.879320v1). RAVA takes FASTQ files (for every sample in your analysis), a metadata sheet (providing info on what day or passage each sample was collected), and a reference genome (your own FASTA and GeneBank files). RAVA outputs an interactive graph (viewable in a web browser), a SNV table for all analyzed samples, and intermediate analysis files.

FASTQ files must be trimmed before running the pipeline.

## Input files
We used only R1 reads for this study. Reads are available at NCBI BioProject PRJNA1306170.

QC filter and trim raw reads with fastp using the following parameters:


```bash
  fastp -i <SAMPLE>.fastq.gz \
    -o "trimmed/<SAMPLE>.fastq.gz" \
    --cut_mean_quality 20 \
    --cut_front --cut_tail \
    --length_required 20 \
    --html "fastp_reports/html/<SAMPLE>_report.html" \
    -j "fastp_reports/json/<SAMPLE>_report.json" \
    --low_complexity_filter \
    --trim_poly_g \
    --trim_poly_x
```


## Running RAVA 

Reference files and an example samplesheet are available in the `references/` and `samplesheets/` directories, respectively.

```bash
nextflow run greninger-lab/RAVA_Pipeline -r v1.0 \
--OUTDIR results/ \
--GFF ./references/CHIKV_181_25_RAVA_ref.gb \
--FASTA ./references/CHIKV_181_25_RAVA_ref.fasta \
--METADATA ./samplesheets/RAVA_samplesheet.csv \
-with-docker ubuntu:18.04 \
-profile Cloud
```

