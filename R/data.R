#' Upper Respiratory Tract samples
#'
#' Human microbiome samples from the upper respiratory tract (URT),
#' taken from a paper by De Boeck et al. It contains nose as well as nasopharynx samples. 
#' Most samples were taken using a swab method, but a minority was taking with the aspirate method.
#'
#' @format ## `urt`
#' A tidytacos object with 1,947 taxa and 207 samples:
#' \describe{
#'   \item{participant}{individual participant ID}
#'   \item{location}{nose (N) or nasopharynx (NF)}
#'   \item{method}{swab (S) or aspirate (A)}
#'   \item{plate}{which of the three plates sequenced the sample is located on}
#'   ...
#' }
#' @source <https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02372/full>
"urt"

#' Phylosphere samples
#'
#' Plant microbiome samples from the phylosphere, taken from a paper by Smets et al.
#'
#' @format ## `leaf`
#' A tidytacos object with 6,044 taxa and 33 samples:
#' \describe{
#'   \item{Plant}{plant species (Blank/Mustard/Wheat)}
#'   \item{Plot}{location of the plant or blank}
#'   \item{Day}{day of sampling}
#'   \item{added_spike_copies}{amount of spike in sequences added}
#'   ...
#' }
#' @source <https://www.sciencedirect.com/science/article/pii/S003807171600050X>
"leaf"