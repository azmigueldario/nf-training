#!/usr/bin/env nextflow

params.greeting = "Hello World!"

include { SPLITLETTERS } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'

include {  START_RUN } from './modules.nf' addParams(foo: 'Óla')
    // <- define new parameter when importing

workflow my_pipeline {
    
    log.info """
        \
        WORKFLOW.NF in Nextflow $nextflow.version
        =========================================
        RunName         : $workflow.runName
        Launch Dir      : $workflow.launchDir
        Work Dir        : $workflow.workDir
        Start time      : $workflow.start
        """
        .stripIndent()
            // add metadata from nextflow and workflow

    // run imported function
    START_RUN()

    // creates a named pipeline
    SPLITLETTERS(Channel.of(params.greeting))
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view{ it }
}

workflow take_pipeline {
    take:
    greeting
        // <- declares one or more input channels

    main:
        // <- if take/emit is used, pipeline is declared under 'main'
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.upper.view {it}
}

workflow emit_pipeline {
    take:
    greeting

    main:
        // <- if take/emit are used, pipeline is declared under 'main'
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())

    emit:
    my_data = CONVERTTOUPPER.out.upper
        // assigns identifier name to outputs
}

workflow unnamed_emit_pipeline {
    take:
    greeting

    main:
    SPLITLETTERS(greeting)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())

    emit:
    CONVERTTOUPPER.out.upper // <- provides output as emit_pipeline.out
}


workflow{
    my_pipeline()
        // A named workflow can be invoked in a latter workflow definition
    take_pipeline(Channel.of(params.greeting))
        // In named workflows, processes can be repeated
        // Using take, we can specify the input(s) while calling the workflow
    emit_pipeline(Channel.of(params.greeting))
    emit_pipeline.out.my_data.view()
        // calling out named ('emit') output
    unnamed_emit_pipeline(Channel.of(params.greeting))
    unnamed_emit_pipeline.out.view()
}