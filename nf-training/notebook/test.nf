sch_test = channel.of(1,2,3)
ch_num = channel.of(1)
//ch_num = channel.value(1)

process SUM {
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
}

// Parsing a JSON

/*
import groovy.json.JsonSlurper

def f = file('../data/meta/regions.json')
def records = new JsonSlurper().parse(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}
*/

