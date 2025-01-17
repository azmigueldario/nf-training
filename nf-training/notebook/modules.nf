/*  Modules are used in a similar way to python in Nextflow's DSL2 (Domain Specific Langugage).
    You can import processes created in a different script to built your main pipeline. Relative paths are usually
    specified to keep it portable. 

    You can specify a new name for the process, or import it twice if 
    necessary (alias required for latter)

        include {PROCESS as ALIAS}  from './module.nf'
    
    For instance, to import multiple modules from the 'hello.nf' script, we could use:

        include {SPLITLETTERS  } from './hello.nf'
        include {CONVERTTOUPPER; REVERSEORDER} from './hello.nf'

*/


process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout emit:upper
        // the 'emit' statement defines name to use in channel scope
    
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

process REVERSEORDER {
    debug true 

    input:
    path x

    output:
    stdout

    script:
    """
    rev $x
    """
}


/*  We can also define custom functions (parameters for Groovy) in our module script.
    This ones can be easily invoked in another script as if they were processes. */

params.foo = "Hola"
params.bar = "Mundo!"

def START_RUN() {
    println "$params.foo $params.bar"
}
