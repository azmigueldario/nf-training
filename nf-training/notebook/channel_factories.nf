// channel.of() factory produces a queue channel

channel
    // 1..23 means all digits starting at 1 and ending at 23
    .of(1..23, 'X', 'Y')
    .view()

// channel.fromList() emits the elements of a list as arguments

list = ['Hello', "World"]
channel
    .fromList(list)
    .view()

// channel.fromPath() emits file(s) matching a glob pattern 

channel
    .fromPath('../data/meta/*.csv')
    .view()

channel
    .fromPath('../data/ggal/**.fq', hidden:true)
    .view()

// channel.fromFilePairs() emits paris of files matching a pattern.
// Channel contains tuples, with a first element ID key and a matching pair of values. 

channel
    .fromFilePairs('../data/ggal/*_{1,2}.fq', hidden:true)
    // must contain at least a star(*) wildcard
    .view()

channel
    .fromFilePairs('../data/ggal/*_{1,2}.fq', flat:true)
    // flat option removes the grouping ID key
    .view()

/*
// channel.fromSRA() looks for matching fastq files in the NCBI SRA repository. 
// Requires an API key to obtain data. File pairs ara managed the same way as if using
// the .fromFilePairs() method.

params.ncbi_api_key = 'Insert your key here'
params.accessions = ['ERR908507', 'ERR908506']

process FASTQC {
    input:
    tuple val(sample_id), path(reads_file)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir -p fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}

workflow {
    reads = channel.fromSRA(params.accessions, apiKey: params.ncbi_api_key)
    fastqc (reads)
}
*/

// .splitText() operator handles multi-line text file parsing into single lines into 
// a channel.

channel
    .fromPath('../data/meta/random.txt')
    .splitText().
    view()

channel
    .fromPath('../data/meta/random.txt')
    .splitText( by:2)
    // .subscribe allows execution of functions after every new value emission 
    .subscribe {
        print it;
        print "----end of the chunk----\n"
    }

channel
    .fromPath('../data/meta/random.txt')
    // closure to transform emited values into upper case font
    .splitText( by:10) { it.toUpperCase() }
    .view()

    // we can also count how many lines
count = 0

channel
    .fromPath('../data/meta/random.txt')
    .splitText()
    .view { "${count++}: ${it.toUpperCase().trim()}" }

  def f = file('../data/meta/random.txt')
  def lines = f.splitText()
  def count=0
  for( String row : lines ) {
    log.info "${count++} ${row.toUpperCase()}"
  }

// .splitCsv() operator is similar to splitText() but for comma separated value files