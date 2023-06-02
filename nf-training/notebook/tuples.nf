/*  Tuples are non-modifiable groups of variables. 
    They can also be used as inputs for nextflow, specifying each of the 
    qualifier and names
    */

reads_ch = Channel.fromFilePairs('../data/ggal/*_{1,2}.fq')

process FOO {
    tag "processing $sample_id"
    debug true 

    input:
    tuple val(sample_id), path(reads)

    output:
    // outputs tuple of a value sample id and a file sample_bam
    tuple val(sample_id), path('sample.bam')

    script:
    """
    echo your_command_here -prefix $sample_id --r1 ${reads[0]} --r2 ${reads[1]} > sample.bam
    """
}


// Exercise naming sample_bam according to sample_id

process EXER1 {
    tag "processing $sample_id"
    debug true 

    input:
    tuple val(sample_id), path(reads)

    output:
    // outputs tuple of a value sample id and a file sample_bam
    tuple val(sample_id), path("${sample_id}.bam")
        // substitution in output assignment requires double quotes

    script:
    """
    echo your_command_here -prefix $sample_id --r1 ${reads[0]} --r2 ${reads[1]} > ${sample_id}.bam
    """
}

workflow {
    // bam_ch = FOO(reads_ch)
    // bam_ch.view { it[1].text }\
        // exercise results
    exer_results = EXER1(reads_ch)
    exer_results.view { "Name: ${it[1].name} -- Content: ${it[1].text}" }
        // To call operators inside double quotes, we use curly braces {}
}