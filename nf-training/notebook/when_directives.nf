/*  The 'when' declaration defines a condition that must be verified
    in order to run the script declaration. 
    It is somehow an 'if statement' and it must turn out in a boolean
    false or true result. 

    To evaluate a regular expression use '=~'
*/

params.dbtype = 'nr'
params.prot = "$projectDir/../data/prots/*.tfa"
proteins = Channel.fromPath(params.prot)

process FIND {
    tag "Processing ${fasta.baseName} file"
    debug true

    input:
    path fasta
    val type

    output:
    path "${fasta.baseName}.txt"

    when:
    fasta.name =~ /^BB11.*/ && type == 'nr'
        // evaluates if fasta.name contains string
        // && (and) evaluates if type is 'nr' 

    script:
    """
    echo blastp -query $fasta -db ${type} > ${fasta.baseName}.txt
    """
}
        /*
        workflow {
            result = FIND(proteins, params.dbtype)
            result.view {it.text}
        }
        */

/*  A 'directive' declaration defines optional settings that 
    determine the execution of the process. They do not affect the 
    semantics of the task itself (input/output/script)

    'directive' declaration is done before 'input' and 'output' declaration.
    'tag' is a common directive that we have used so far. 

    They are commonly used to assign the computing resources or
    to configure other parameters. 
*/

// 'publishDir' directive specifies where to store files for long-term use. 
// By default, all files are kept inside the 'work/' folder. 

params.outdir = 'my_results'
params.prot2 = "$projectDir/../data/prots/*.tfa"
proteins = Channel.fromPath(params.prot2)

process BLASTSEQ{
    publishDir "$params.outdir/bam_files", mode: 'copy'
        // mode: copy leaves a backup inside the 'work/' folder

    input:
    path fasta

    output:
    path ('*.txt')

    script:
    """
    echo blastp -query $fasta -db nr > ${fasta}_result.txt
    """
}

// To manage subdirectories, we can specify matching patterns.

params.reads = "$projectDir/data/reads/*_{1,2}.fq.gz"
params.outdir = 'my_results'
samples_ch = Channel.fromFilePairs(params.reads, flat:true)

process FOO{
    publishDir "$params.outdir/$sampleId/", pattern: '*.fq'
    publishDir "$params.outdir/$sampleId/counts", pattern: '*_counts.txt'
    publishDir "$params.outdir/$sampleId/outlooks", pattern: '*_outlooks.txt'

    input:
    tuple val(sampleid), path('sample1.fq.gz'), path('sample2.fq.gz')

    output:
    path "*"

    when:
    params.outdir != "$projectDir"

    script:
    """
    sample1.fq.gz zcat > sample1.fq
    sample2.fq.gz zcat > sample2.fq

    awk '{S++}END{print s/4}' sample1.fq > sample1_counts.txt
    awk '{S++}END{print s/4}' sample2.fq > sample2_counts.txt

    head -n 50 sample1.fq > sample1_outlook.txt
    head -n 50 sample1.fq > sample1_outlook.txt
    """
    

}


workflow {
    blast_ch = BLASTSEQ(proteins)
    blast_ch.view { "${it.baseName} contains: ${it.text}"}
        // when applies a function/operator, use {}
}