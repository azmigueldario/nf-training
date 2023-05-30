/* 
 * pipeline input parameters 
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir="results"

"""
println "reads: $params.reads"
println "outdir: $params.outdir"  
"""

// this prints the log.info below
log.info """\
        nextflow/rnaseq-nf
        =============================================
        reads                           : ${params.reads}
        transcript                      : ${params.transcript}
        multiqc                         : ${params.multiqc}
        outdir                          : ${params.outdir}
        """
        .stripIndent()