process CreateGff_Genbank_RAVA { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3
	
    input:
      val(GENBANK)
	  file PULL_ENTREZ
	  file WRITE_GFF

    output: 
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    script:
    """
    #!/bin/bash
	python3 ${PULL_ENTREZ} ${GENBANK}
	/usr/local/miniconda/bin/bwa index lava_ref.fasta
	python3 ${WRITE_GFF}

    """
}

// Pulls --GENBANK fasta
process Pull_Genbank { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1
	
    input:
      val(GENBANK)
      file CONTROL_FASTQ
	  file PULL_ENTREZ

    output: 
      file "lava_ref.fasta"
	  file "CONTROL.fastq"

    script:
    """
    #!/bin/bash
    
	set -e 
	echo ${CONTROL_FASTQ}
    
	python3 ${PULL_ENTREZ} ${GENBANK}
	# Avoiding filename collision during run_pipeline process 
	mv ${CONTROL_FASTQ} CONTROL.fastq
    """
}

process CreateGff_RAVA { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3
	
    input:
      val(GENBANK)
	  file PULL_ENTREZ
	  file WRITE_GFF
	  file FASTA 
	  file GFF

    output: 
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    script:
    """
    #!/bin/bash

	grep -v "mature_peptide" ${GFF} > lava_ref.gff
	grep "mature_peptide" ${GFF} | sed "s/,mature_peptide//g" > mat_peptides.txt
	mv ${FASTA} lava_ref.fasta
	#mv ${GFF} lava_ref.gff
	#Creates empty txt file
	touch ribosomal_start.txt
	#touch mat_peptides.txt
	cp lava_ref.fasta consensus.fasta
	/usr/local/miniconda/bin/bwa index lava_ref.fasta

    """
}

process CreateGFF { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1
	
    input:
      val(GENBANK)
      file CONTROL_FASTQ
	  file PULL_ENTREZ
	  file WRITE_GFF
	  file FASTA
	  file GFF

    output: 
      file "lava_ref.fasta"
      file "lava_ref.gff"
	  file "CONTROL.fastq"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    script:
    """
    #!/bin/bash
    
	set -e 
	echo ${CONTROL_FASTQ}
    
	grep -v "mature_peptide" ${GFF} > lava_ref.gff
	grep "mature_peptide" ${GFF} | sed "s/,mature_peptide//g" > mat_peptides.txt
	mv ${FASTA} lava_ref.fasta
	touch ribosomal_start.txt
	
	# Avoiding filename collision during run_pipeline process 
	 mv ${CONTROL_FASTQ} CONTROL.fastq
    """
}

process Align_first_sample {
	container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    errorStrategy 'retry'
    maxRetries 3

    input:
	file CONTROL_FASTQ
	file "lava_ref.fasta"

	output: 
	file "CONTROL.bam"

	shell:
	'''
	#!/bin/bash

	# Align each sample to consensus fasta.
	/usr/local/bin/bbmap.sh in=!{CONTROL_FASTQ} outm=CONTROL.bam ref=lava_ref.fasta local=true -Xmx6g
	'''
}

process Process_first_sample {
	container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
	file "CONTROL.bam"
	file VCFUTILS
	file "lava_ref.fasta"

	output:
	tuple file("CONTROL.vcf.gz"), file("CONTROL.vcf.gz.tbi")
	file "lava_ref.fasta.fai"

	shell:
	'''
	#!/bin/bash
	/usr/local/miniconda/bin/samtools faidx lava_ref.fasta
	/usr/local/miniconda/bin/samtools sort -@ !{task.cpus} CONTROL.bam > CONTROL.sorted.bam
    /usr/local/miniconda/bin/samtools index CONTROL.sorted.bam

	# Unwrap fasta
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' lava_ref.fasta > a.tmp && mv a.tmp lava_ref.fasta

	splitnum=$(($(($(tail -n +2 lava_ref.fasta | awk '{print length}')/!{task.cpus}))+1))
	# Generates pileup that VCF can be called off of later.
	perl !{VCFUTILS} splitchr -l $splitnum lava_ref.fasta.fai | \\
            xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                "/usr/local/miniconda/bin/bcftools mpileup \\
                    -f lava_ref.fasta -r {} \\
                    --count-orphans \\
                    --no-BAQ \\
                    --max-depth 50000 \\
                    --max-idepth 500000 \\
                    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                CONTROL.sorted.bam | /usr/local/miniconda/bin/bcftools call -A -m -Oz - > tmp.{}.vcf.gz"
        
        # Concatenate parallelized vcfs back together
        gunzip tmp*vcf.gz
        mv tmp.*\\:1-* CONTROL_catted.vcf
        for file in tmp*.vcf; do grep -v "#" $file >> CONTROL_catted.vcf; done

        cat CONTROL_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | /usr/local/miniconda/bin/bcftools norm -m -any > CONTROL_pre_bcftools.vcf
		/usr/local/miniconda/bin/bcftools filter -i 'IMF > 0.5 || (DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' --threads !{task.cpus} CONTROL_pre_bcftools.vcf -o CONTROL_pre2.vcf
        /usr/local/miniconda/bin/bcftools norm --check-ref s --fasta-ref lava_ref.fasta -Ov CONTROL_pre2.vcf > CONTROL_pre3.vcf

		# pull out header
        grep "#" CONTROL_pre3.vcf > CONTROL.vcf
        # get rid of long sgRNAs that are called due to fake depths
        grep -v "#" CONTROL_pre3.vcf | awk -F'\t' 'length($4) <60 { print }' >> CONTROL.vcf

        # Index and generate consensus from vcf with majority variants
        /usr/local/miniconda/bin/bgzip CONTROL.vcf
        /usr/local/miniconda/bin/tabix CONTROL.vcf.gz

	'''
}

// Generate final consensus from first sample vcf.
process Generate_consensus {
    container "quay.io/biocontainers/bcftools:1.14--hde04aa1_1"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple file("CONTROL.vcf.gz"), file("CONTROL.vcf.gz.tbi")
		file "lava_ref.fasta"
		file "lava_ref.fasta.fai"
    output:
        file "consensus.fasta"

    script:
    """
        cat lava_ref.fasta | /usr/local/bin/bcftools consensus CONTROL.vcf.gz > consensus.fasta
    """
}

// Uses accession number specified by --GENBANK to create our own GFF
process CreateGFF_Genbank { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1
	
    input:
      file "consensus.fasta"
	  file WRITE_GFF

    output: 
      file "lava_ref.gff"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    script:
    """
    #!/bin/bash
    
	set -e 
	python3 ${WRITE_GFF}
    """
}


// Prep for downstream steps.
// Indexes and prepares consensus fasta, and generates prerequisite files for Annovar.
process Alignment_prep { 
    container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
	  file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    output: 
	file "consensus.fasta"
	file "AT_refGene.txt"
	file "AT_refGeneMrna.fa"
    file "lava_ref.gff"
	file "ribosomal_start.txt"
	file "mat_peptides.txt"
	file "AT_refGeneMrna.fa"

	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: 'AT_refGeneMrna*'

    // Code to be executed inside the task
    script:
    """
    #!/bin/bash

	# Unwrap fasta
	awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' consensus.fasta > a.tmp && mv a.tmp consensus.fasta
	
	# Preparatory steps for Annovar downstream
	gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes
	retrieve_seq_from_fasta.pl --format refGene --seqfile consensus.fasta AT_refGene.txt --out AT_refGeneMrna.fa 
    """
}

// Aligns all samples to consensus fasta.
process Align_samples { 
	container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    errorStrategy 'retry'
    maxRetries 3

    input:
	tuple file(R1), val(PASSAGE)
	file "consensus.fasta"

	output: 
	tuple file(R1), file("*.bam"), val(PASSAGE)

	shell:
	'''
	#!/bin/bash
	echo aligning "!{R1}"

	# Align each sample to consensus fasta.
	/usr/local/bin/bbmap.sh in=!{R1} outm=!{R1}.bam ref=consensus.fasta local=true -Xmx6g
	'''
}

// Generates genomecov files and pileups.
process Process_samples { 
	container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
	tuple file(R1), file("${R1}.bam"), val(PASSAGE)
	file VCFUTILS
	file ATREF
	file ATREF_MRNA
	file "consensus.fasta"

	output: 
	tuple file(R1), file("${R1}.bam"), val(PASSAGE)
	file "${R1}.genomecov"
	file "*exonic_variant_function" optional true
	tuple file(R1), file("*.bam"), file( "*.exonic_variant_function.samp"), val(PASSAGE)
	file "${R1}.vcf"

	publishDir "${params.OUTDIR}vcf_files", mode: 'copy', pattern: '*.vcf'

	shell:
	'''
	#!/bin/bash
	# Creates genomecov file from BAM so we can generate coverage graphs later.
	echo sample\tposition\tcov > !{R1}.genomecov
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  !{R1}.bam >> !{R1}.genomecov

	/usr/local/miniconda/bin/samtools faidx consensus.fasta
	/usr/local/miniconda/bin/samtools sort -@ !{task.cpus} !{R1}.bam > !{R1}.sorted.bam
    /usr/local/miniconda/bin/samtools index !{R1}.sorted.bam

	splitnum=$(($(($(tail -n +2 consensus.fasta | awk '{print length}')/!{task.cpus}))+1))
	# Generates pileup that VCF can be called off of later.
	perl !{VCFUTILS} splitchr -l $splitnum consensus.fasta.fai | \\
            xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                "/usr/local/miniconda/bin/bcftools mpileup \\
                    -f consensus.fasta -r {} \\
                    --count-orphans \\
                    --no-BAQ \\
                    --max-depth 50000 \\
                    --max-idepth 500000 \\
                    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                !{R1}.sorted.bam | /usr/local/miniconda/bin/bcftools call -A -m -Oz - > tmp.{}.vcf.gz"
        
        # Concatenate parallelized vcfs back together
        gunzip tmp*vcf.gz
        mv tmp.*\\:1-* !{R1}_catted.vcf
        for file in tmp*.vcf; do grep -v "#" $file >> !{R1}_catted.vcf; done

        cat !{R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | /usr/local/miniconda/bin/bcftools norm -m -any > !{R1}.vcf
		awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","1/0",$10)gsub("1/1","1/0",$11)}1\' !{R1}.vcf > !{R1}_p.vcf

		file="!{R1}""_p.vcf"
		convert2annovar.pl -withfreq -format vcf4 -includeinfo !{R1}_p.vcf > !{R1}.avinput 
		annotate_variation.pl -outfile !{R1} -v -buildver AT !{R1}.avinput .
		mv !{R1}.exonic_variant_function !{R1}.exonic_variant_function.samp	
	'''
}

// Initializes proteins.csv - list of protein names and locations - from our generated GFF.
process Pipeline_prep { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file blank_ignore
		file "lava_ref.gff"
		file "consensus.fasta"
		file INITIALIZE_PROTEINS_CSV

	output: 
		file 'merged.csv'
		file 'proteins.csv'
	
	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: 'proteins.csv'

	script:
	"""
	#!/bin/bash
	# Creates header for final csv.
	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv
	# Creates list of protein names and locations (proteins.csv) based on GFF annotations.
	python3 ${INITIALIZE_PROTEINS_CSV}
	"""
}

// Extract variants for all samples.
process Extract_variants { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(R1), file(BAM), file(EXONICVARIANTS), val(PASSAGE)
		file METADATA
		file AT_REFGENE_MRNA
		file RIBOSOMAL_LOCATION
		file PROTEINS_CSV
		file FIX_VARIANTS_FILE
		file MAT_PEPTIDE_LOCATIONS
		file MAT_PEPTIDE_ADDITION

	output:
		tuple file("${R1}.csv"), val(PASSAGE), file("reads.csv"), file(R1) optional true
		tuple file(R1), val(PASSAGE) optional true
		file R1
		file "${R1}.reads.csv"
		file "${R1}.csv"
		file "${R1}.genomecov"

	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: '*.genomecov'
	
	shell:
	'''
	#!/bin/bash
	echo !{R1}
	# Creates genomecov files for genome coverage graphs later.
	echo 'sample	position	cov' > !{R1}.genomecov 
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam !{BAM} >> !{R1}.genomecov
	# reads.csv from all processes will be merged together at end 
	printf !{R1}"," > reads.csv
	/usr/local/miniconda/bin/samtools flagstat !{BAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv
	#awk -F":" '($26+0)>=1{print}' !{EXONICVARIANTS}> !{R1}.txt
	grep -v "UNKNOWN" !{EXONICVARIANTS} > variants.txt

    # Sorts by beginning of mat peptide
    sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
    # Adds mature peptide differences from protein start.
    python3 !{MAT_PEPTIDE_ADDITION}
    
	python3 !{FIX_VARIANTS_FILE} -name !{R1} -fasta !{AT_REFGENE_MRNA} -metadata !{METADATA}
	
	cp reads.csv !{R1}.reads.csv
	'''
}

// Generates LAVA visualization plots for whole genome and for each protein across samples.
process Generate_output { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file R1 
		file READS_CSV
		file SAMPLE_CSV
		file MERGED_CSV
		file PROTEINS_CSV
		file GENOMECOV
		file VCF		
		file GENOME_PROTEIN_PLOTS
		file PALETTE
		file MAT_PEPTIDE_LOCATIONS
		file MAT_PEPTIDE_ADDITION
		file METADATA
		val(TITLE)

	output:
		file "*.html"
		file "*.csv"
	
	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: '*.py'
	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: 'reads.csv'
	publishDir "${params.OUTDIR}", mode: 'copy', pattern: 'LAVA_plots.html'
	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: 'visualization.csv'
	publishDir "${params.OUTDIR}visualization", mode: 'copy', pattern: 'mat_peptides_add*'

	script:

	"""
	#!/bin/bash
	ls -lah

	# Sorts by beginning of mat peptide
    sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
    # Adds mature peptide differences from protein start.
    python3 ${MAT_PEPTIDE_ADDITION}

	cat *.reads.csv > reads.csv
	echo "Sample,AminoAcidChange,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage,MatPeptide" > visualization.csv
	
	# Order file by what is in metadata file
	for i in \$(tail -n +2 ${METADATA}); do j=\$(echo \$i | rev | cut -d'/' -f1 | rev | cut -d',' -f1); tail -n +2 \$j'.csv' >> visualization.csv; done

	python3 ${GENOME_PROTEIN_PLOTS} visualization.csv proteins.csv reads.csv . "${TITLE}"
	
	"""
} 
