// Define parameters

params.ref_genome = './data/genome.fa'
params.known_variants = './data/known_variants.vcf.gz'
params.blacklist = './data/blacklist.bed'
params.outdir = "$launchDir/results"

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

workflow{
    genome_ch = channel.fromPath(params.ref_genome)
    //genome_dir = file(params.ref_genome).getParent()
    index_ch = SAMTOOLS_INDEX(genome_ch)
    dictionary_ch = PICARD_DICTIONARY(genome_ch)
    star_index_ch = STAR_INDEX(genome_ch)
}


workflow.onComplete {
    println "Done, workflow finished with status: ${ workflow.success ? 'Completed' : 'Failed' }"
}

/*
workflow.onComplete {
    println "Created the following documents: ${index_ch.index.view()}, ${dictionary_ch.dictionary.view()}"
}
process TEST {
    debug true
    tag "Processing $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    stdout

    script:
    """
    echo $sample_id --readFilesIn ${reads[0]} ${reads[1]}
    """
}

workflow {
    TEST(reads_ch).view()

}
*/