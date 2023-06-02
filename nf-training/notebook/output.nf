/*  The output declaration block defines the channels to send the process results.
    You can only have one output block, but multiple channels in it. */

// 'value' qualifier defines a value in the script context. 

methods = ['prot', 'dna', 'rna']

process FOO_0 {
    debug true

    input:
    val x

    output:
    val x //  <- must match a value in the script

    script:
    """
    echo $x > file
    """
}

// 'path' qualifier specifies one or more files produced by the process

process RANDOMNUM {
    debug true
    
    output:
    path 'result.txt'

    script: // <- single quotes, so $ calls a nextflow variable
    '''
    # RANDOM is nextflow variable for random number
    echo $RANDOM > result.txt
    '''

}

        /*
        workflow {
                // val output
            receiver_ch = FOO_0(Channel.of(methods))
            receiver_ch.view { "Received: $it"}
                // path output
            receiver2_ch = RANDOMNUM()
            receiver2_ch.view { "Received: " + it.text }
        }
        */

// Specifying options for path qualifier output

process EX1 {
    output:
        // applies the hidden option to output file
    path 'result.txt', hidden: true // <- requires colon ':'
        

    '''
    echo 'another new line' >> result.txt
    '''
}

process EX2 {
    
    output:
    // parenthesis are mandatory () when multiple outputs qualifiers are assigned
    // and you need to specify the options for one or several of them
    tuple path('last_result.txt'), path('result.txt', hidden: true)

    '''
    echo 'another new line' >> result.txt
    echo 'another new line' > last_result.txt
    '''
}

// multiple output files and wildcards

process SPLITLETTERS {
    input:
    val x

    output:
    path 'file_*'

    """
        # '$x' in single quotes is expanded completely
    printf '$x' | split -b 1 - file_
    """
            // -b means byte_count
            // - <file_prefix>
}
        /*
        workflow{
            input_ch = SPLITLETTERS (channel.of("Hola mundo bonito!"))
            input_ch
                .flatMap()
                // flatMap transforms list of files into channel with individual files
                .view { "file: ${it.name} ==> ${it.text}"}
        }
        */

// Dynamic file names as output

species = ['Pseudomonas_aeruginosa', 'Mycobacterium_tuberculosis', 'Klebsiella_pneumoniae']
sequences = ['AATGGGAT', 'CCTGGGAC', 'CCCTCCAGG']

Channel.fromList(species)
    .set { species_ch }

process ALIGN {
    tag "Processing $organism"
    debug true

    input:
    val organism
    val seq

    output:
    path "${organism}.txt"

    script:
    """
    echo align -in $seq > ${organism}.txt
    """
}

workflow {
    genomes = ALIGN(species_ch, sequences)
    genomes.view { "Contents: ${it.text}" }
}