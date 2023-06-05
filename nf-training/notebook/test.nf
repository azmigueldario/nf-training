/*
        ch_test = channel.of(1,2,3)
        ch_num = channel.of(1)
        //ch_num = channel.value(1)

        process SUM {
            debug true

            input:
            val x
            val y

            output:
            stdout

            script:
            """
            echo \$(($x + $y))
            """
        }

        // The .first() operator transforms the queue channel into a value channel 

        workflow {
            SUM(ch_test, ch_num.first()).view()
            // ch_test.view()
            println "$workflow.projectDir"
        }
*/

// Parsing a JSON

/*
        import groovy.json.JsonSlurper

        def f = file('../data/meta/regions.json')
        def records = new JsonSlurper().parse(f)


        for( def entry : records ) {
        log.info "$entry.patient_id -- $entry.feature"
        }
*/

/*
params.dbtype = 'nr'
params.prot = '../data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

params.reads = "$projectDir/../data/reads/*_{1,2}.fq.gz"


process FIND {
  debug true

  input:
  path fasta
  val type

  when:
  fasta.name =~ '/^BB11.{*}/ && type == 'nr'

  script:
  """
  echo blastp -query $fasta -db nr
  """
}
*/
    /*
    workflow {
    result = FIND(proteins, params.dbtype)
    Channel.fromFilePairs(params.reads, flat:true)
        .view()
    }
    */

params.greeting = "Hello world! My name is "
params.name = "Miguel"
message = (params.greeting + params.name)
println(message + "ssss")
//message = params.greeting.concat(params.name)
greetin_ch = Channel.of(params.greeting + params.name)
greetin_ch.view()




