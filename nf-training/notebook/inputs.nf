/*  The `val` qualifier receives data of any type as input, 
    it is accessed in the script with the input name */

num = channel.of(1, 2, 3)

process BASICEXAMPLE{
    debug true

    input:
    val x

    script:
    // input is called by name, here is `x`
    """
    echo process job $x
    """
}

/*  The `path` qualifier  allows the handling of file values.
    It acepts one or more files, and the use of wildcard symbols.
    We can access it on the script using the assigned `name` or a variable reference. */

reads = Channel.fromPath ('../data/ggal/*.fq')
reads2 = channel.fromPath('../data/reads/*.fq.gz')
exercise_ch = channel.fromPath('../data/ggal/*_1.fq')


process FOO {
    debug true // <- prints to std output for debugging 

    input:
    path 'sample.fastq' // <- with name
    path sample // <- variable reference

    script:
    """
    ls sample.fastq
    ls $sample
    """
}

process FOO2 {
    debug true 

    input:
    path sample

    script:
        // expects multiple files 
    """
    ls -lha $sample
    """

}

process EXER1 {
    tag "Concat all files" // <- debug true would not work as no stdout is produced

    input:
    path exer_file

    output:
    path 'top_10_lines'

    script:
    """
    cat $exer_file > merged_file
    head -n 20 merged_file > top_10_lines
    """
}

    /*
    workflow{
        // myrun = BASICEXAMPLE(num)
        // result = FOO(reads, reads2)
        // FOO2(reads2.collect()) // <- collect merges all files from path input
        concat_ch = EXER1(exercise_ch.collect())
        concat_ch.view()
    }
    */

/*  Combine input channels. When one of the channels is empty, the process is 
    terminated, even if a process has more values. */

ch1 = channel.of(1, 2, 3)
ch2 = channel.of('a', 'b', 'c')
    // channels with different number of elements
input1 = channel.of(1, 2)
input2 = channel.of('a', 'b', 'c', 'd', 'e')


process FOO3 {
    debug true 

    input:
    val x
    val y

    script:
        // x and y will be assigned in the order they are specified in workflow
    """
    echo $x and $y
    """
}

process DIFFERENT_CARDINALITY {
    debug true 

    input:
    val 'first_input'
    val 'second_input'

    script:
    """
    echo cardinality: ${first_input} and ${second_input}
    """
}

// Channel.value() can be recicled several times as in the following process

input3 = channel.value(90)

process CHANNEL_VALUE_RECYCLE {
    debug true
    
    input:
    val x
    val y

    script:
    """
    echo recycle: $x and $y
    """
}

    /*
    workflow {
        FOO3(ch1, ch2) // <- ch1 assigned to $x and ch2 to $y
        DIFFERENT_CARDINALITY(input1, input2) // <- different number of values
        CHANNEL_VALUE_RECYCLE(input3, input2)
    }
    */

// Exercise to create a process that repeats a value channel

transcriptome_file = "$baseDir/../data/ggal/transcriptome.fa"
readsEx_ch = channel.fromPath('../data/ggal/*_1.fq')

process EXER2 {
    debug true

    input:
    path transcriptome
    path read_1

    output:
    path result

    script:
    """
    echo salmon -q -r1 $read_1 -t $transcriptome > result
    """
}

workflow {
    concat_ch = EXER2(transcriptome_file, readsEx_ch)
    concat_ch.view()
}

