# Efficacy of the 4’-Fluorouridine nucleoside analog against Nipah virus in the Syrian hamster model
This repository contains the code and data used in the publication "Efficacy of the 4’-Fluorouridine nucleoside analog against Nipah virus in the Syrian hamster model" by Escaffre et al. (2025).

## What is RAVA?
RAVA is derived from LAVA, but for non-longitudinal sequence data. RAVA takes FASTQ files (for every sample in your analysis), a metadata sheet (providing info on what day or passage each sample was collected), and a reference genome either by accession number or your own FASTA and GFF. RAVA outputs an interactive graph (viewable in a web browser), a mutation table for all analyzed samples, and intermediate analysis files.

FASTQ files must be trimmed before running the pipeline.

## Input files
All raw read files for this study can be found in BioProject PRJNA1214874.
QC filter and trim raw reads with fastp using the following parameters:

```
  fastp -i <raw_reads_file_name>.fastq.gz \
    -o "trimmed/<reads_trimmed_file_name>.fastq.gz" \
    --cut_mean_quality 20 \
    --cut_front --cut_tail \
    --length_required 20 \
    --html "fastp_reports/html/<file_name>_report.html" \
    -j "fastp_reports/json/<file_name>_report.json" \
    --low_complexity_filter \
    --trim_poly_g \
    --trim_poly_x
```

You can subsample afg46 (P0) reads using the following command:
```bash
seqtk sample -s112 <(gunzip -c afg46.fastq) 2000000 > afg46_2Msubsample.fastq
```

## Running RAVA on P0 stock (sample afg46) used for animal infection

```bash
nextflow run greninger-lab/RAVA_Pipeline -r publications_2025_01_01 \
--OUTDIR nf-out/P0_sample-afg46_results/ \
--GFF ./references/rNiVB-Gluc-P2A-eGFP-AG_RAVA_ref.gb \
--FASTA ./references/rNiVB-Gluc-P2A-eGFP-AG_RAVA_ref.fasta \
--METADATA ./samplesheets/<your_samplesheet>.csv \
-with-docker ubuntu:18.04 \
-profile Cloud
```

Where `<your_samplesheet>.csv` is a CSV file with the following columns:
```csv
Sample,Passage
./path/to/afg46_2Msubsample.fastq.gz,afg46
```

## Running RAVA on lung and brain samples

```bash
nextflow run greninger-lab/RAVA_Pipeline -r publications_2025_01_01 \
--OUTDIR nf-out/lung_and_brain_results/ \
--GFF ./references/rNiVB-Gluc-P2A-eGFP-afg46-AG_RAVA_ref.gb \
--FASTA ./references/rNiVB-Gluc-P2A-eGFP-afg46-AG_RAVA_ref.fasta \
--METADATA ./samplesheets/<your_samplesheet>.csv \
-with-docker ubuntu:18.04 \
-profile Cloud
```

Where `<your_samplesheet>.csv` is a CSV file with the following columns (include only samples that you are interested in - one sample per line):

```csv
Sample,Passage
./path/to/gr2a_bl_lung.fastq.gz,gr2a_bl_lung
./path/to/gr2a_wh_lung.fastq.gz,gr2a_wh_lung
./path/to/gr2a_yl_lung.fastq.gz,gr2a_yl_lung
...
```

- See example samplesheet in `samplesheets/` directory.
- All reference files used in the analysis can be found in the `references/` directory.
- Pre-generated RAVA plots for afg46 (P0), lung and brain samples can be found in the `RAVA_plots/` directory.
- Full lists of SNV results generated for afg46 (P0), lung and brain samples can be found in the `SNV_results/` directory.