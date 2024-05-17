// Uses accession number specified by --GENBANK to create our own GFF (lava_ref.gff) for a consensus fasta
// generated from the alignment of "Passage 0" sample to reference fasta.

process CreateGFF {

    container "quay.io/vpeddu/lava_image:latest"

    input:
		file PULL_ENTREZ
		file WRITE_GFF
		file FASTA
		file GFF

    output:
		file "*.fasta"
		file "*.gff"
		file "ribosomal_start.txt"
		file "mat_peptides.txt"

    script:

    """

    #!/bin/bash

	sed '/^>/! s/U/T/g' ${FASTA} > lava_ref_T.fasta

	if [[ ${GFF} == *.gff ]]

		then

			grep -v "mature_peptide" ${GFF} > lava_ref.gff
			grep "mature_peptide" ${GFF} | sed "s/,mature_peptide//g" > mat_peptides.txt
			#Creates empty txt file
			touch ribosomal_start.txt
			/usr/local/miniconda/bin/bwa index ${FASTA}

		else

			/usr/local/miniconda/bin/bwa index ${FASTA}
			mv ${GFF} lava_ref.gbk
			python3 ${WRITE_GFF} ${FASTA}
	
	fi

	FASTA_BASENAME=\$(basename ${FASTA} .fasta)

	# Convert U to T fasta 
	sed '/^>/! s/U/T/g' ${FASTA} > lava_ref_T.fasta
	mv lava_ref_T.fasta \${FASTA_BASENAME}_RAVA.fasta
	mv lava_ref.gff \${FASTA_BASENAME}_RAVA.gff			 

    """

}

// Aligns all samples to consensus fasta and removes duplicates if --DEDUPLICATE specified.
// Also generates genomecov files and pileups.

process Align_samples {

	container "quay.io/vpeddu/lava_image:latest"

    input:
	tuple file(R1), val(PASSAGE)
	file FASTA
	val DEDUPLICATE

	output:
	tuple file(R1), file("*.pileup"), file("*.bam"), val(PASSAGE)
	file "${R1}.genomecov"
	file "${R1}.bam"

	shell:

	'''

	#!/bin/bash

	/usr/local/miniconda/bin/bwa index !{FASTA}
	/usr/local/miniconda/bin/samtools faidx !{FASTA}
	gatk CreateSequenceDictionary -R !{FASTA} --VERBOSITY ERROR --QUIET true

	echo aligning "!{R1}"

	# Align each sample to consensus fasta.
	/usr/local/miniconda/bin/bwa mem -t !{task.cpus} -M -R \'@RG\\tID:group1\\tSM:!{R1}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -L [2,2] -B 6 !{FASTA} !{R1} > !{R1}.sam

	# Sorts SAM.
	java -jar /usr/bin/picard.jar SortSam INPUT=!{R1}.sam OUTPUT=!{R1}.bam SORT_ORDER=coordinate VERBOSITY=ERROR

	# Removes duplicates (e.g. from library construction using PCR) if --DEDUPLICATE flag specified.
	if !{DEDUPLICATE}
		then
			echo "Deduplicating !{R1}"
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${R1}_dedup.bam > ${R1}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=!{R1}.bam VERBOSITY=ERROR

	# Creates genomecov file from BAM so we can generate coverage graphs later.
	echo sample\tposition\tcov > !{R1}.genomecov
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  !{R1}.bam >> !{R1}.genomecov

	# Generates pileup that VCF can be called off of later.
	/usr/local/miniconda/bin/samtools mpileup -B --max-depth 500000 -f !{FASTA} !{R1}.bam > !{R1}.pileup

	'''

}

// Initializes proteins.csv - list of protein names and locations - from our generated GFF.
process Pipeline_prep {

    errorStrategy 'retry'
    maxRetries 3
	container "quay.io/vpeddu/lava_image:latest"

	input:
		file GFF
		file INITIALIZE_PROTEINS_CSV

	output:
		file 'merged.csv'
		file 'proteins.csv'

	script:

	"""

	#!/bin/bash

	# Creates header for final csv.
	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage,Reverse_Complement" > merged.csv

	# Creates list of protein names and locations (proteins.csv) based on GFF annotations.
	python3 ${INITIALIZE_PROTEINS_CSV} ${GFF}

	"""

}

// Generates VCF for all the samples and converts to .avinput for Annovar.

process Create_VCF {

	container "quay.io/vpeddu/lava_image:latest"

	input:
		tuple file(R1), file(R1_PILEUP), file(BAM), val(PASSAGE)
		file ANNOCAR
		file FASTA
		file GFF 

	output:
		file "*exonic_variant_function" optional true
		tuple file(R1), file("*.bam"), file( "*.exonic_variant_function.samp"), val(PASSAGE)
		file "${R1}.vcf"

	shell:

	'''

	#!/bin/bash

	ls -latr
	echo Analyzing variants in sample !{R1}

	# here for file passthrough (input -> output)
	mv !{BAM} !{BAM}.bam

	# Generates VCF outputting all bases with a min coverage of 2.
	cat !{R1_PILEUP} | java -jar /usr/local/bin/VarScan mpileup2cns --validation 1 --output-vcf 1 --min-coverage 2 --min-var-freq 0.001 --p-value 0.99 --min-reads2 1 --strand-filter 1 > !{R1}.vcf

	# Fixes ploidy issues.
	awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{R1}.vcf > !{R1}_p.vcf

	# Converts VCF to .avinput for Annovar.
	file="!{R1}""_p.vcf"
	convert2annovar.pl -withfreq -format vcf4old -includeinfo !{R1}_p.vcf > !{R1}.avinput

	python3 annoCAR.py !{R1}.avinput !{GFF} !{FASTA}

	mv !{R1}.exonic_variant_function !{R1}.exonic_variant_function.samp

	'''

}

// Extract variants for all samples.

process Extract_variants {

    errorStrategy 'retry'
    maxRetries 3
	container "quay.io/vpeddu/lava_image:latest"

	input:
		tuple file(R1), file(BAM), file(EXONICVARIANTS), val(PASSAGE)
		file METADATA

	output:
		tuple file("${R1}.csv"), val(PASSAGE), file("reads.csv"), file(R1) optional true
		tuple file(R1), val(PASSAGE) optional true

	shell:

	'''

	#!/bin/bash

	# Creates genomecov files for genome coverage graphs later.
	echo 'sample	position	cov' > !{R1}.genomecov
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam !{BAM} >> !{R1}.genomecov

	# reads.csv from all processes will be merged together at end
	printf !{R1}"," > reads.csv
	/usr/local/miniconda/bin/samtools flagstat !{BAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

	awk -F":" '($26+0)>=1{print}' !{EXONICVARIANTS}> !{R1}.txt

	grep "SNV" !{R1}.txt > a.tmp
	grep "stop" !{R1}.txt >> a.tmp
	mv a.tmp !{R1}.txt
	SAMPLE="$(awk -F"," -v name=!{R1} '$1==name {print $2}' !{METADATA})"
	echo $SAMPLE

	awk -v name=!{R1} -v sample=!{PASSAGE} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$44+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$41","sample","$(NF)}' !{R1}.txt > !{R1}.csv

	'''

}

// Checks for multi-nucleotide mutations and prints out warning message.
// Currently LAVA does not handle complex mutations and instead annotates it as such for manual review.

process Annotate_complex {

    errorStrategy 'retry'
    maxRetries 3
	container "quay.io/vpeddu/lava_image:latest"

	input:
		tuple file(SAMPLE_CSV), val(PASSAGE), file("reads.csv"), file(R1)
		file ANNOTATE_COMPLEX_MUTATIONS

	output:
		file R1
		file "${R1}.complex.log"
		file "${R1}.reads.csv"
		file SAMPLE_CSV

	script:

	"""

	#!/bin/bash

	# Checks for complex mutations and prints a warning message.
	python3 ${ANNOTATE_COMPLEX_MUTATIONS} ${SAMPLE_CSV} ${PASSAGE}

	# Renaming files to avoid file collision
	mv complex.log ${R1}.complex.log
	mv reads.csv ${R1}.reads.csv

	"""

}

// Generates LAVA visualization plots for whole genome and for each protein across samples.

process Generate_output {

	container "quay.io/vpeddu/lava_image:latest"

	input:
		file R1
		file COMPLEX_LOG
		file READS_CSV
		file SAMPLE_CSV
		file MERGED_CSV
		file PROTEINS_CSV
		file GENOMECOV
		file VCF
		file RIBOSOMAL_LOCATION
		file MAT_PEPTIDE_LOCATIONS
		file MAT_PEPTIDE_ADDITION
		file RIBOSOMAL_SLIPPAGE
		file GENOME_PROTEIN_PLOTS
		file PALETTE
		file BAM
		val(NAME)
		file FASTA

	output:
		file "*.html"
		file "*.log"
		file "visualization.csv"
		file "reads.csv"
		file "vcf_files"
		file "genomecov"
		file "all_files"
		file "bam_files"
		publishDir params.OUTDIR, mode: 'copy'

	script:

	"""

	#!/bin/bash

	ls -lah
	cat merged.csv > final.csv

	# Takes fastq.gz and fastq
	if ls *.gz &>/dev/null
	then
		cat *.fastq.gz.csv >> final.csv
	else
		cat *.fastq.csv >> final.csv
	fi
	
	# Gets rid of non-SNPs
	grep -v "transcript" final.csv > a.tmp && mv a.tmp final.csv
	grep -v "delins" final.csv > a.tmp && mv a.tmp final.csv

	# Sorts by beginning of mat peptide
	sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt

	# Adds mature peptide differences from protein start.
	python3 ${MAT_PEPTIDE_ADDITION}
	rm mat_peptides.txt

	# Corrects for ribosomal slippage.
	python3 ${RIBOSOMAL_SLIPPAGE} final.csv proteins.csv
	awk NF final.csv > a.tmp && mv a.tmp final.csv
	cat *.reads.csv > reads.csv
	cat *.log > complex.log

	if [ "${NAME}" == "false" ]; then
		HTML_NAME=\$(basename "${FASTA}" .fasta)
	else
		HTML_NAME="${NAME}"
	fi

	# TODO error handling @ line 669-683 of lava.py
	python3 ${GENOME_PROTEIN_PLOTS} visualization.csv proteins.csv reads.csv . "Plot" -name \${HTML_NAME} 
	mkdir vcf_files
	mv *.vcf vcf_files

	mkdir genomecov
	mv *.genomecov genomecov

	mkdir bam_files
	mv *.bam bam_files 

	mkdir all_files

	cp -r *.txt all_files

	"""

}