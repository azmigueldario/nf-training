// Define parameters

params.ref_genome = './data/genome.fa'
params.known_variants = './data/known_variants.vcf.gz'
params.blacklist = './data/blacklist.bed'
params.outdir = "$launchDir/results"
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar"

// Create channels from parameters

known_variants_ch = Channel.fromPath('./data/known_variants.vcf.gz')
reads_ch = Channel.fromFilePairs('./data/reads/*_{1,2}.fastq.gz')

process SAMTOOLS_INDEX {
    debug true
    tag "processing $genome"

    input:
    path genome

    output: 
    path "${genome}.fai", emit:index

    script:
    """
    samtools faidx $genome
    """
}

process PICARD_DICTIONARY {
    debug true
    tag "Creating genome dictionary"
    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    input:
    path genome

    output:
    path "${genome.baseName}.dict", emit: dictionary 

    script:
    """
    picard CreateSequenceDictionary R= $genome O= "${genome.baseName}.dict"
    """
}

process STAR_INDEX {
    debug true
    tag "Creating STAR index for $genome"
    publishDir "${params.outdir}"

    input:
    path genome

    output:
    path 'genome_dir'

    script:
    """
    mkdir genome_dir
    # must be inside work dir

    STAR --runMode genomeGenerate \
        --genomeDir genome_dir \
        --genomeFastaFiles ${genome} \
        --runThreadN ${task.cpus}
    """
}

process FILTER_BLACKLISTED{
    debug true

    input:
    path known_variants
    path blacklisted

    output:
    tuple path("filtered_variants.vcf.gz"), path("filtered_variants.vcf.gz.tbi")
    
    script:
    """
    vcftools    --gzvcf $known_variants \
                -c                                      `# option '-c' is necessary in newer versions of vcftools` \
                --exclude-bed $blacklisted \
                --recode | bgzip -c \
                > filtered_variants.vcf.gz
    
    tabix filtered_variants.vcf.gz
    """
}

process STAR_ALIGNMENT {

    input:
    path genome
    path genome_index
    tuple val(meta_id), path(reads)

    output:
    tuple val(meta_id), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')

    script:
    """
    # align reads to genome
    STAR --genomeDir $genome_index \
         --readFilesIn $reads \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999

    # second pass (improve alignments using splice junctions table 
    # and create a new index)
    mkdir genomeDir
    STAR --runMode genomeGenerate \
         --genomeDir genomeDir \
         --genomeFastaFiles $genome \
         --sjdbFileChrStartEnd SJ.out.tab \
         --sjdbOverhang 75 \
         --runThreadN ${task.cpus}

    # Final read alignments
    STAR --genomeDir genomeDir \
         --readFilesIn $reads \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:$meta_id LB:library PL:illumina PU:machine SM:GM12878

    # Index the BAM file
    samtools index Aligned.sortedByCoord.out.bam
    """
}

process SPLIT_N {
    tag "Processing $meta_id"
    debug true
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"
        /*  
            1.  Evaluates if containerEngine is 'singularity' and if the 'task.ext.singularity_pull_docker_container'
                option was not specified.
            2.  If the conditions apply, then it downloads singularity image 'https://depot.galaxy..'
            3.  If not, it downloads from the docker repository 'quay.io/biocontainers...'  
        */

    input:
    path genome_ref
    path genome_dict
    path genome_index
        // <- Index and dictionary are necessary for 'SplitNCigarReads' but not specified in command
    tuple val(meta_id), path(bam), path(bai)

    output:
    tuple val(meta_id), path ('split_bam')

    script:
    """
    gatk SplitNCigarReads \
        --reference $genome_ref \
        --input $bam \
        --output split_bam \
        --refactor-cigar-string
    """
}

process GATK_RECALIBRATE {
    debug true
    tag "Running $meta_id"
    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0'

    input:
    path ref_genome
    path genome_dict
    path genome_index
    tuple path(known_variants), path(known_variants_index)
    tuple val(meta_id), path (split_bam)

    output:
    tuple val(meta_id), path('final_bam')


    script:
    """
        # run recalibrator
    gatk BaseRecalibrator \
        --known-sites $known_variants \
        --reference $ref_genome \
        --input $split_bam \
        --output rnaseq_grp

    gatk ApplyBQSR \
          -R $ref_genome \
          -I $split_bam \
          --bqsr-recal-file rnaseq_grp \
          -O final_bam
    """

}

process REMOVE_MULTIALIGN {
    input:
    tuple val(meta_id), path (final_bam)

    output:
    tuple val(meta_id), path ('final_uniq_bam'), path('final_uniq_bam.bai')

    script:
    """
        # Select only unique alignments, no multimaps
    (samtools view -H final_bam; samtools view final_bam | grep -w 'NH:i:1') \
    |samtools view -Sb - > final_uniq_bam
  
        # Index the BAM file
    samtools index final_uniq_bam
    """
}

workflow{
    // define input channels
    genome_ch = channel.fromPath(params.ref_genome)
    blacklisted_ch = channel.fromPath(params.blacklist)
    known_variants_ch = channel.fromPath(params.known_variants)

    // preparation of inputs
    samtools_index_ch = SAMTOOLS_INDEX(genome_ch)
    dictionary_ch = PICARD_DICTIONARY(genome_ch)
    star_index_ch = STAR_INDEX(genome_ch)
    filtered_variants_ch = FILTER_BLACKLISTED(known_variants_ch, blacklisted_ch)

    
    final_alignment_ch = STAR_ALIGNMENT(genome_ch, star_index_ch, reads_ch)
    split_bam_ch=SPLIT_N(genome_ch, dictionary_ch, samtools_index_ch, final_alignment_ch)
    GATK_RECALIBRATE(genome_ch, dictionary_ch, samtools_index_ch, filtered_variants_ch, split_bam_ch)
    final_uniq_bam_ch = REMOVE_MULTIALIGN(GATK_RECALIBRATE.out)
}



workflow.onComplete {
    println "Done!. Workflow finished with status: ${ workflow.success ? 'Completed' : 'Failed' }"
}
