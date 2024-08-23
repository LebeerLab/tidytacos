#' Construct a phylogenetic tree from the sequences in the taxa table of a tidytacos object.
#'
#' `add_tree()` performs a multiple sequence alignment on the ASV sequences in
#' the taxa table and then constructs a phylogenetic tree.
#' The tree is added to the tidytacos object under a variable called `tree`.
#'
#' @param ta A tidytacos object.
#' @param sequence_var The name of the column in the taxa table that contains the sequence
#' @param outgroup_taxon_id The taxon_id of the outgroup to root the tree. If NULL, the tree will be unrooted.
#' @param aln_args An optional list of arguments to pass to the [DECIPHER::AlignSeqs()] function.
#' @param tree_fit_args An optional list of arguments to pass to the [phangorn::optim.pml()] function.
#' @return A tidytacos object with a phylotree slot.
#' @examples 
#' # filter taxa to speed up calculation time
#' \dontrun{urt_sub <- urt %>% 
#'   filter_taxa(taxon_id %in% head(urt$taxa$taxon_id))
#' urt_sub_tree <- add_tree(urt_sub, 
#'                          sequence_var="sequence", 
#'                          outgroup_taxon_id="t4", 
#'                          aln_args=list(verbose=FALSE)
#'                         )
#' # our new rooted tree
#' urt_sub_tree$tree
#' 
#' urt_sub_unrooted <- add_tree(urt_sub,aln_args=list(verbose=FALSE))
#' # an unrooted version of the tree
#' urt_sub_unrooted$tree}
#' @importFrom methods as
#' @importFrom stats update
#' @export
add_tree <- function(ta, sequence_var=sequence, outgroup_taxon_id=NULL, aln_args=list(),tree_fit_args=list()) {
    force_optional_dependency('DECIPHER', instructions = 'It can be installed with:\n\tBiocManager::install("DECIPHER")')
    force_optional_dependency('phangorn', instructions = 'It can be installed with:\n\tinstall.packages("phangorn")')
    force_optional_dependency('ape', instructions = 'It can be installed with:\n\tinstall.packages("ape")')

    sequence <- rlang::enquo(sequence_var)
    
    # Get the sequences
    if (!rlang::as_name(sequence) %in% names(ta$taxa)) {
        stop(paste("No sequences found under the column", sequence_var, "in the taxa table."))
    }
    seqs <- ta$taxa %>%
      pull(!!sequence) %>%
      Biostrings::DNAStringSet()
    names(seqs) <- ta$taxa$taxon_id
    
    # Align them
    #aln <- DECIPHER::AlignSeqs(seqs,anchor=NA) %>% 
    aln_args["anchor"] <- NA
    aln <- eval(rlang::call2(
        "AlignSeqs",seqs, !!!aln_args, .ns="DECIPHER"
        )) %>%
      as("matrix")
    
    ph_aln <- phangorn::phyDat(aln, "DNA")
    dist_m <- phangorn::dist.ml(ph_aln)
    treeNJ <- phangorn::NJ(dist_m)
    
    fit <- phangorn::pml(treeNJ, data=ph_aln)
    fitGTR <- update(fit, k=4, inv=0.2)

    if (!any(c("optInv", "optGamma", "rearrangements") %in% names(tree_fit_args))) {
        # use default values
        tree_fit_args["optInv"] <- TRUE
        tree_fit_args["optGamma"] <- TRUE
        tree_fit_args["rearrangements"] <- "stochastic"
    }
    fitGTR <- eval(rlang::call2("optim.pml",
        fitGTR, control = phangorn::pml.control(trace=0), !!!tree_fit_args, .ns="phangorn"
        ))
                 
    ta$tree <- phyloseq::phy_tree(fitGTR$tree)
    if (!is.null(outgroup_taxon_id)) {
        ta$tree <- ape::root(ta$tree, outgroup_taxon_id, resolve.root=TRUE)
    }
    ta
}

# internal function
ensure_tree <- function(ta, rooted=FALSE){
    force_optional_dependency('ape', instructions = 'It can be installed with:\n\tinstall.packages("ape")')
    if(rooted && (!"tree" %in% names(ta)|| !ape::is.rooted(ta$tree))){
        stop(paste("No rooted tree found in the tidytacos object.",
        "\nPlease run add_tree(outgroup_taxon_id='tx') first, where tx",
        "is the taxon_id of the outgroup."))
    } else if (!"tree" %in% names(ta)) {
        stop("No tree found in the tidytacos object. Please run add_tree() first.")
    }
}
