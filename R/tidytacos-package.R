#' @keywords internal
"_PACKAGE"

#' tidytacos: Functions to manipulate and visualize amplicon count data.
#'
#' tidytacos is an R package for the analysis of amplicon count data:
#' abundances of amplicon sequences (either clustered in OTUs or exact variants)
#' in samples. The package builds on the tidyverse created by Hadley Wickham:
#' the data are stored in tidy tables where each row is an observation and each
#' column a variable. In addition, the package supplies a set of "verbs":
#' functions that take a tidytacos object as first argument and also return
#' a tidytacos object. Not all functionality is currently implemented in the
#' form of verbs, but this will soon be remediated.
#'
#' @docType package
#' @name tidytacos
#'
#' @section Author(s):
#' Stijn Wittouck \email{wittouck_stijn@@hotmail.com}
#'

## usethis namespace: start
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @import vegan
#' @importFrom tibble tibble
#' @importFrom stats dist
globalVariables(c(
    "sample_id", "samples",
    "taxon_id", "taxa", "taxon", "taxon_name", "taxon_name_color",  
    "count", "rel_abundance", "presence", "alpha_metrics",
    "sample_clustered", "prevalence", "cluster", "total_absolute_abundance",
    "spike_count", "total_count", "mean_rel_abundance", "spike_abundance",
    "total_count", "prevalence", "total_density",
    ".", "packageVersion", "packageDescription", "ord1", "ord2",
    "colorRampPalette"
))
## usethis namespace: end

# Class -----------------------------------------------------------------------
#' An S4 class to store a grouped tidytaco with \link{group_samples}
#'
#' This S4 class stores the output from \link{group_samples}.
#'
#' @slot groups the different grouped tidytacos.
#' @name grouped_taco
#' @importFrom methods new
setClass(
    "grouped_taco",
    slots=c(
        tacos="list",
        group_cols="character"
    )
)
