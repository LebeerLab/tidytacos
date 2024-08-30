#' Construct a phylogeny from ASV sequences
#'
#' `add_tree()` performs a multiple sequence alignment on the ASV sequences in
#' the taxa table and then constructs a phylogenetic tree.
#' The tree is added to the tidytacos object under a variable called `tree`.
#'
#' @param ta A tidytacos object.
#' @param sequence_var The name of the column in the taxa table that contains 
#'   the ASV sequences
#' @param aln_args An optional list of arguments to pass to the 
#'   [DECIPHER::AlignSeqs()] function
#' @param tree_fit_args An optional list of arguments to pass to the 
#'   [phangorn::optim.pml()] function
#' @return A tidytacos object with a phylotree slot
#' 
#' @examples
#' \dontrun{
#' # filter taxa to speed up calculation time
#' urt_sub <- urt %>% filter_taxa(taxon_id %in% head(urt$taxa$taxon_id))
#' # infer tree
#' urt_sub_tree <- add_tree(urt_sub, 
#'                          sequence_var="sequence", 
#'                          aln_args=list(verbose=FALSE)
#'                         )
#' # inspect tree
#' urt_sub_tree$tree
#' }
#' 
#' @importFrom methods as
#' @importFrom stats update
#' @export
add_tree <- 
    function(ta, sequence_var=sequence, aln_args=list(), tree_fit_args=list()) {
        
    force_optional_dependency('DECIPHER', 
      instructions = 'It can be installed with:\n\tBiocManager::install("DECIPHER")')
    force_optional_dependency('phangorn', 
      instructions = 'It can be installed with:\n\tinstall.packages("phangorn")')
    force_optional_dependency('ape', 
      instructions = 'It can be installed with:\n\tinstall.packages("ape")')

    sequence <- rlang::enquo(sequence_var)
    
    # Get the sequences
    if (!rlang::as_name(sequence) %in% names(ta$taxa)) {
        stop(paste(
          "No sequences found under the column", sequence_var, 
          "in the taxa table."
        ))
    }
    seqs <- ta$taxa %>%
      pull(!!sequence) %>%
      Biostrings::DNAStringSet()
    names(seqs) <- ta$taxa$taxon_id
    
    # Align sequences
    aln_args["anchor"] <- NA
    aln <- 
      eval(rlang::call2(
        "AlignSeqs", seqs, !!!aln_args, .ns="DECIPHER"
        )) %>%
      as("matrix")
    
    # Infer NJ tree
    ph_aln <- phangorn::phyDat(aln, "DNA")
    dist_m <- phangorn::dist.ml(ph_aln)
    treeNJ <- phangorn::NJ(dist_m)
    
    # Optimize branch lengths using GTR model
    fit <- phangorn::pml(treeNJ, data=ph_aln)
    fitGTR <- update(fit, k=4, inv=0.2)
    settings <- c("optInv", "optGamma", "rearrangements")
    defaults <- list(TRUE, TRUE, "stochastic")
    for (i in 1:length(settings)) {
        if (! settings[i] %in% tree_fit_args) {
            tree_fit_args[settings[i]] <- defaults[[i]]
        }
    }
    fitGTR <- eval(rlang::call2("optim.pml",
        fitGTR, control = phangorn::pml.control(trace=0), !!!tree_fit_args, 
        .ns="phangorn"
        ))  
    ta$tree <- phyloseq::phy_tree(fitGTR$tree)
    
    # Root tree
    ta$tree <- phangorn::midpoint(ta$tree)
        
    ta
        
}

# internal function
ensure_tree <- function(ta, rooted=FALSE){
    force_optional_dependency('ape', 
      instructions = 'It can be installed with:\n\tinstall.packages("ape")')
    if(rooted && (!"tree" %in% names(ta)|| !ape::is.rooted(ta$tree))){
        stop(paste("No rooted tree found in the tidytacos object.",
        "\nPlease run add_tree(outgroup_taxon_id='tx') first, where tx",
        "is the taxon_id of the outgroup."))
    } else if (!"tree" %in% names(ta)) {
        stop("No tree found in the tidytacos object. Please run add_tree() first.")
    }
}
