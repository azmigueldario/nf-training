#!usr/bin/env nextflow

// Using modules, this script is much more straightforward to understand. 

params.greeting = 'Hello world! My name is '
params.name = "Miguel"
//greeting_ch = Channel.of(params.greeting + params.name)
    // <- The '+' operator concatenates strings without space


include { SPLITLETTERS as DIVIDELETTERS   } from './modules.nf'
include { CONVERTTOUPPER as CHANGETOMAYUS } from './modules.nf'
    // <- use an alias for the process in this script, similar to python
include { SPLITLETTERS; CONVERTTOUPPER } from './modules.nf'
    // <- load multiple processes from the same script, separate with ';'

    /*
    workflow {
        letters_ch1 = SPLITLETTERS(greeting_ch)
        results_ch1 = CONVERTTOUPPER(letters_ch1.flatten())
        results_ch1.view()

        letters_ch2 = DIVIDELETTERS(greeting_ch).flatten()
        letters_ch2 = CHANGETOMAYUS(letters_ch2)
        letters_ch2.view()
    }
    */

// Output definition in workflow

workflow  {
    greeting_ch = Channel.of(params.greeting + params.name)
        // eliminate channel definitions, a channel is implied every line
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out[0].view()
        // when multiple output is available, each can be accessed with [ ]
    CONVERTTOUPPER.out.upper.view { it }
        // can also used named output when identified with 'emit'

// Below I use pipes to send output from one process to the other

    Channel.of (params.greeting) |
        DIVIDELETTERS |
        flatten() |
        CHANGETOMAYUS |
        view
}