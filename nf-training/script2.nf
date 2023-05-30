/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that creates a binary index 
 * given the transcriptome file
 */
process INDEX {
    cpus 4
    input:
    path transcriptome
     
    output:
    path 'salmon_index' 

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

workflow {
    index_ch = INDEX(params.transcriptome)
    index_ch.view()
}