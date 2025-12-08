# Impact of Mutations Affecting 4’-Fluorouridine Susceptibility on Fitness and Treatment Outcomes for Venezuelan Equine Encephalitis Virus

This repository contains the code and data used in the publication **"Impact of Mutations Affecting 4’-Fluorouridine Susceptibility on Fitness and Treatment Outcomes for Venezuelan Equine Encephalitis Virus"** by **Wong** et al. (2026).

## What is RAVA?

RAVA is derived from [LAVA](https://www.biorxiv.org/content/10.1101/2019.12.17.879320v1). RAVA takes FASTQ files (for every sample in your analysis), a metadata sheet (providing info on what day or passage each sample was collected), and a reference genome (your own FASTA and GeneBank files). RAVA outputs an interactive graph (viewable in a web browser), a SNV table for all analyzed samples, and intermediate analysis files.

FASTQ files must be trimmed before running the pipeline.

## Input files

We used only R1 reads for this study. Reads are available at NCBI BioProject PRJNA1306184.

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

The consensus sequence of VEEV-EGFP P0 (`240213-GSU-RC-025-VEEV-EGFP_Consensus.gb`) was generated in Geneious Prime by aligning reads from the P0 input sample to the plasmid reference sequence (`plasmid_map_240213-GSU-RC-025-VEEV-EGFP.gb`) using majority voting. P0 consensus sequence was then used as a reference in RAVA analysis. Reference files and an example samplesheet are available in the `references/` and `samplesheets/` directories, respectively.

```bash
nextflow run greninger-lab/RAVA_Pipeline -r v1.0 \
--OUTDIR results/ \
--GFF ./references/240213-GSU-RC-025-VEEV-EGFP_Consensus.gb \
--FASTA ./references/240213-GSU-RC-025-VEEV-EGFP_Consensus.fasta \
--METADATA ./samplesheets/RAVA_samplesheet.csv \
-with-docker ubuntu:18.04 \
-profile Cloud
```

## Strand analaysis

To investigate if specific mutations were located on the same viral genome or on different genomes, a strand analysis was performed using [Positional Nucleotide Profiler v0.2](https://github.com/DariiaVyshenska/positional_nuc_profiler) using indexed BAM files that were generated during RAVA analysis.

Example command:

```bash
positional_nuc_profiler my_file_sorted_indexed.bam . 6139 6261 6274
```

## Notes

Folder `selected_ooutput_files` contains full lists of detected SNVs across all samples as well as SNV metadata.