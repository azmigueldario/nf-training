// splitCsv() is an operator similar to splitText() but for comma
// separated values. 

channel
    .fromPath('../data/meta/patients_1.csv')
    .splitCsv()
    .view  { row -> "${row[0]}, ${row[3]}"} // <- row is a list object

    // if file has header, it can be specified and used to call values by column
channel
    .fromPath('../data/meta/patients_1.csv')
    .splitCsv(header:true)
    .view  { row -> "${row.patient_id}, ${row.num_samples}"} // <- row is a list object

    // can also provide custom header names
channel
    .fromPath('../data/meta/patients_1.csv')
    .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'])
    .view  { row -> "${row.col1}, ${row.col4}, ${row.col5}"} 

    // multiple .csv files can be processed simultaneously
channel
    .fromPath('../data/meta/patients_*.csv') // <- matches a glob pattern
    .splitCsv(header: true)
    .view  { row -> "${row.patient_id}\t${row.num_samples}"} // <- \t is a tab separator

    // can be used without channel definition
def f = file('../data/meta/patients_1.csv') // <- expects a single file
    def lines = f.splitCsv()
    for ( List row : lines) {
        log.info "${row[0]} -- ${row[2]}"
    }

// tsv files can be parsed using .splitCsv() too with the option sep:'\t'
channel
    .fromPath('../data/meta/regions.tsv', checkIfExists:true)
    // use 'sep' as an option to parse TAB separated files
    .splitCsv(sep:'\t', skip:1)
    .view { row -> "$row[0]" }