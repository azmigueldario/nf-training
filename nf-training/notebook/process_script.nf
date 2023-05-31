params.data = "World"

process EXAMPLE {

  script:
  """
  echo 'Hello world!\nHola mundo!\nCiao mondo!\nHallo Welt!' > file
  cat file | head -n 1 | head -c 5 > chunk_1.txt
  gzip -c chunk_1.txt  > chunk_archive.gz
  """
}

process PYSTUFF {

    // specify output to standard output
    output:
    stdout

    // by default script is in bash, but can use other languages
    script:
    """
    #!/usr/bin/env python

    x = 'Hello'
    y = 'world!'
    print ("%s - %s" % (x,y))
    """
}

process FOO {
    
    script:
    """
    echo Hello params.data
    """
}

process FOO2 {

    script:
    // escape an environment variable when using double quotes
    """
    echo "The current working directory is \$PWD"
    """
}

process BAR {
    script:
    """
    echo $PATH | tr : '\\n'
    """
}

workflow {
    PYSTUFF() | view {"$it"} // <- output can be seen with view operator
    EXAMPLE()
    FOO()
    FOO2()
    BAR()
}