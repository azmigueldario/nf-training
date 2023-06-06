#!usr/bin/env nextflow

params.greeting = 'Hello world! My name is '
params.name = "Miguel"
greetin_ch = Channel.of(params.greeting + params.name)
    // <- The '+' operator concatenates strings without space

process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout
    
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

process REVERSEORDER {
    debug true 

    input:
    path x

    output:
    stdout

    script:
    """
    rev $x
    """
}

workflow {
    letters_ch1 = SPLITLETTERS(greetin_ch)
    results_ch1 = CONVERTTOUPPER(letters_ch1.flatten())
    results_ch1.view()
}