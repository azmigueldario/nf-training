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