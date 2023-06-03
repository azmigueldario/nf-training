/*  Operators are something like methods in python. 
    They can be used to transform a channel's emission, or apply rules. 

*/

//  randomSample() takes a sample of the specified 'n'

Channel
    .of(1..1000)
    .randomSample(12, 439) // <- second value is seed for reproducibility
    //.view()

//  unique() removes duplicates from the channel when emitted

Channel
    .of (1, 2, 3, 6, 6, 6, 7, 7, 8, 1, 4)
    .unique { it % 4 } // <- defines characteristics of uniqueness, picks five only
    //.unique()
    //.view()

//  collect() operator obtains all items emitted by a channel into a list and returns the 
//  list as a single emission. 

Channel
    .fromPath("$projectDir/../data/ggal/*.fq")
    .collect()
    //.view()

Channel
    .of('hello', 'oi', 'hola')
    .collect { it.length() } // <- returns length of each emission
    //.view()

Channel
    .of(['hello', 'hola'], ['oi', 'bonjour', 'hi'], ['konnichiwa', 'namaste', 'ciao'] )
    .collect(flat: true) // <- puts everything in the same list
    //.view()

/*  concat() concatenates the items emitted by a channel in the same order
    they are in the operator arguments. */

a_ch = Channel.of('a', 'b', 'c')
b_ch = Channel.of(1, 2, 3)
c_ch = Channel.of('p', 'q')

    //c_ch.concat(b_ch, a_ch).view() // <- emits c completely, then appends c, and then a_ch

/*  join() creates a new channel that joins the contents of two channels
    according to a matching key. By default it is the first element in each
    emitted item. */

left_ch = channel.of (['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right_ch = Channel.of(['Z', 6], ['Y', 5], ['X', 2])

    /*
    left_ch.join(right_ch).view()
        // merges based on first element emmited each time, in this case the letter
    left_ch.join(right_ch, by: 1, remainder: true).view() 
        // <- merge by item[1] or second one
        // <- remainder emits unmatched pairs
    */

// The mix() operator combines the items emitted by two or more channels into a single one

c1 = Channel.of( 1, 2, 3 )
c2 = Channel.of( 'a', 'b' )
c3 = Channel.of( 'z' )

    // c1.mix(c2,c3)
        // .subscribe onNext: { println it }, onComplete: { println 'Done' }
        // <- emits in any order
        // <- subscribe specifies a function after every emission

//  The map() operator applies a function of your choosing to every emitted item.
//  The function is expressed with a closure inside curly brackets {}

Channel
    .of(1, 2, 3, 4, 5)
    .map {it * it}
    .subscribe onNext: {println it}, onComplete: { println "Done." }

Channel
    .fromPath('../data/ggal/*_1.fq')
        // can be used to emit tuples too
    .map { pair -> [pair.getBaseName(), pair] }
    .view { name, pair -> "> $name: $pair"}