#!/usr/bin/env nextflow
/*
========================================================================================
                        LAVA
========================================================================================
Longitudinal Analysis of Viral Alleles
#### Homepage / Documentation
https://github.com/greninger-lab/RAVA_Pipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    RAVA: Reference-based Analysis of Viral Alleles
    Usage:
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/lava \\
        --METADATA      Required argument: A two column csv - the first column is the
                        path to all the fastqs you wish to include in your analysis.
                        All fastqs that you want to include need to be specified in
                        this file AND be located in the folder from which you are
                        running lava. The second column is the temporal seperation
                        between the samples. This is unitless so you can input
                        passage number, days, or whatever condition your experiment
                        happens to have. [REQUIRED]
        --OUTDIR        Output directory [REQUIRED]
        --FASTA         Specify a reference fasta to map samples to. This must be the 
                        same fasta as your annotation file (.gb or .gff). This option 
                        must be used with the --GFF flag to specify the protein 
                        annotations relative to the start of this fasta. 
        --GFF           Specify a reference .gff, .gb or .gbk file with the protein annotations for
                        the reference fasta supplied with the --FASTA flag. This option
                        must be paired with the --FASTA flag.
        --AF            pecify an allele frequency percentage to cut off
                        - with a minimum of 1 percent - in whole numbers. default = ' 
        --ALLELE_FREQ   Specify an allele frequency percentage to cut off - with a
                        minimum of 1 percent - in whole numbers.
        --DEDUPLICATE   Optional flag, will perform automatic removal of PCR
                        duplicates via DeDup.
        --NAME          Changes graph name
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.OUTDIR= false
params.GFF = 'False'
params.FASTA = 'NO_FILE'
params.DEDUPLICATE = 'false'
params.ALLELE_FREQ = 'NO_VAL'

params.NAME = 'false'

METADATA_FILE = file(params.METADATA)

/*
 * Import the processes used in this workflow
 */

include { CreateGFF } from './Modules.nf'
include { Align_samples } from './Modules.nf'
include { Pipeline_prep } from './Modules.nf'
include { Create_VCF } from './Modules.nf'
include { Extract_variants } from './Modules.nf'
include { Annotate_complex } from './Modules.nf'
include { Generate_output } from './Modules.nf'

// Staging python scripts
WRITE_GFF = file("$workflow.projectDir/bin/write_gff.py")

INITIALIZE_PROTEINS_CSV = file("$workflow.projectDir/bin/initialize_proteins_csv.py")
ANNOTATE_COMPLEX_MUTATIONS = file("$workflow.projectDir/bin/Annotate_complex_mutations.py")
MAT_PEPTIDE_ADDITION = file("$workflow.projectDir/bin/mat_peptide_addition.py")
RIBOSOMAL_SLIPPAGE = file("$workflow.projectDir/bin/ribosomal_slippage.py")
GENOME_PROTEIN_PLOTS = file("$workflow.projectDir/bin/genome_protein_plots.py")
PALETTE = file("$workflow.projectDir/bin/palette.py")
ANNOCAR = file("$workflow.projectDir/bin/annoCAR.py")

// Error handling for input flags
//if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR")
    exit(1)
}

// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
    params.OUTDIR = "${params.OUTDIR}/"
}

input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.Sample), (row.Passage)) }

// Throws exception if paths in METADATA are not valid
Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{row-> (file(row.Sample, checkIfExists:true))}.ifEmpty{error "Check metadata file"}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    
    log.info nfcoreHeader()

        CreateGFF (
            PULL_ENTREZ,
            WRITE_GFF,
            file(params.FASTA),
            file(params.GFF)
        )

        Align_samples (
            input_read_ch,
            CreateGFF.out[0],
            params.DEDUPLICATE
        )

        Pipeline_prep (
            CreateGFF.out[1],
            INITIALIZE_PROTEINS_CSV
        )

        Create_VCF (
            Align_samples.out[0],
            ANNOCAR,
            CreateGFF.out[0],
            CreateGFF.out[1]
        )

        Extract_variants (
            Create_VCF.out[1],
            METADATA_FILE
        )

        Annotate_complex(
            Extract_variants.out[0],
            ANNOTATE_COMPLEX_MUTATIONS
        )

        Generate_output(
            Annotate_complex.out[0].collect(),
            Annotate_complex.out[1].collect(),
            Annotate_complex.out[2].collect(),
            Annotate_complex.out[3].collect(),
            Pipeline_prep.out[0],
            Pipeline_prep.out[1],
            Align_samples.out[1].collect(),
            Create_VCF.out[2].collect(),
            CreateGFF.out[2],
            CreateGFF.out[3],
            MAT_PEPTIDE_ADDITION,
            RIBOSOMAL_SLIPPAGE,
            GENOME_PROTEIN_PLOTS,
            PALETTE, 
            Align_samples.out[2].collect(),
            params.NAME,
            file(params.FASTA)
        )
}

def nfcoreHeader() {

    return """
                       ooO
                     ooOOOo
                   oOOOOOOoooo
                 ooOOOooo  oooo
                /vvv\\
               /V V V\\
              /V  V  V\\
             /         \\            oh wow  look at these alleles
            /           \\          /         /
          /               \\   	  o          o
__       /                 \\     /-   o     /-
/\\     /                     \\  /\\  -/-    /\\
                                    /\\


        `7MM"'"Mq.        db `7MMF'   `7MF' db
          MM   `MM.      ;MM:  `MA     ,V  ;MM:
          MM   ,M9      ,V^MM.  VM:   ,V  ,V^MM.
          MMmmdM9      ,M  `MM   MM.  M' ,M  `MM
          MM  YM.      AbmmmqMA  `MM A'  AbmmmqMA
          MM   `Mb.   A'     VML  :MM;  A'     VML
        .JMML. .JMM..AMA.   .AMMA. VF .AMA.   .AMMA.

                      Version 6.0.0
    """.stripIndent()

}
