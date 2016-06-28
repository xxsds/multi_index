


pretty_input <- function( input_col ) {
    input_col <- basename(as.character(input_col))
    input_col <- gsub(".data$","",input_col)
    input_col <- gsub("Clueweb09-Full.SimHash","\\\\cluesim",input_col)
    input_col <- gsub("Clueweb09-Full.OddSketch","\\\\clueodd",input_col)
    input_col <- gsub("lsh_sift_64.hash","\\\\siftlsh",input_col)
    input_col <- gsub("mlh_sift_64.hash","\\\\siftmlh",input_col)
    return (input_col)
}

pretty_query <- function ( query_col ) {
    query_col <- basename(as.character(query_col))
    query_col <- gsub("\\.query$","",query_col)
    query_col <- gsub("(.*\\.)(.*)","\\2",query_col)
    if ( sum( data$query %in% c("existing","real","random") ) == 0 ){
        query_col <- "query"
    }
    return (query_col)
}
