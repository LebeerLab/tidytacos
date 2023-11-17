#' Determine the correlation between the distance of the counts in a tidytacos object and a sample variable, multiple sample variables or another matrix.
#'
#' This function performs a mantel test using the dissimilarity matrix  
#' of the tidytacos object supplied and a second distance matrix generated from user input.
#'
#' @param ta A tidytacos object.
#' @param comparison A distance to compare against. This can be any of the following:
#'  - The name of the variable in the sample table to use for comparison
#'  - A list of names of variables in the sample table.
#'  - A distance matrix.
#' 
#' @return The mantel test statistics
#'
#' @export
perform_mantel_test <- function(ta, comparison, ...) {

    dmatrix <- ta %>% counts_matrix()
    # Matrix get sorted on sample, so we need to sort the samples in the ta object
    ta$samples <- ta$samples %>% arrange(sample)
    if (length(dmatrix[,1]) < length(ta$samples$sample_id)) {
        warning("Empty samples found, ignoring them in analysis")
        ta <- ta %>% remove_empty_samples()
    }

    if (typeof(comparison) == "double") {
        return(vegan::mantel(dmatrix, comparison, ...))
    } 

    if (length(comparison) > 1) { 
        comparison <- ta$samples %>% select(comparison)
        return(mantel_test_list(dmatrix, comparison, ...))
    }

    comparison <- ta$samples %>% dplyr::pull(comparison)
    return(mantel_test_vector(dmatrix, comparison, ...))
}

mantel_test_vector <- function(dmatrix, vector, ...) {
    
    if (typeof(vector) == "character") {
        vector <- as.numeric(as.factor(vector))
    }
    dist.vector <- dist(vector)
    vegan::mantel(dmatrix, dist.vector, ...)
}

mantel_test_list <- function(dmatrix, parameters, ...) {
    
    force_optional_dependency("fastDummies")
    parameters <- fastDummies::dummy_cols(parameters) %>% 
      select_if(is.numeric)
    parameters <- scale(parameters, center=T, scale=T)
    d.param <- dist(parameters)
    print(typeof(d.param))
    vegan::mantel(dmatrix, d.param, ...)   
}


