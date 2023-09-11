#' Add relative abundance to abundance table
#'
#' \code{add_rel_abundance} adds relative abundance to the abundance table of a
#' tidytacos object.
#'
#' This function adds the relative abundance per sample to the abundance table
#' of a tidytacos object under the variable name "rel_abundance".
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE)
#'
#' # Add relative abundance
#' data <- data %>%
#'  add_rel_abundance()
#'
#' @export
add_rel_abundance <- function(ta) {

  # add relative abundance to abundance table
  ta$counts <- ta$counts %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = count / sum(count)) %>%
    ungroup()

  # return ta object
  ta

}

#' Add absolute abundance to counts table
#'
#' \code{add_absolute_abundance} adds absolute abundance to the counts table of a
#' tidytacos object.
#'
#' This function adds the absolute abundance per taxon and per sample to the counts table
#' of a tidytacos object under the variable name "absolute_abundance". It also adds total absolute abundances to the samples table.
#'
#' @param ta tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the samples table which indicates how much spike was added per sample, e.g. 16S rRNA gene copy numbers added to the DNA extraction tube.
#' @param material_sampled The amount of material from which DNA was extracted, e.g gram of soil. This parameter encourages researchers to consider that absolute abundances are only meaningful if they can be translated into densities.
#' 
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'  c(1500, 1300, 14, 280, 356, 9),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2", "taxon3")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE)
#'
#' # Add total abundance
#' data <- data %>%
#'  add_absolute_abundance(spike_taxon="t3")
#'
#' @export
add_absolute_abundance <- function(ta, spike_taxon, spike_added=1000, material_sampled=1) {
	  
  # if total_counts and relative abundances not present: add temporarily
  total_counts_tmp <- ! "total_counts" %in% names(ta$samples)
  if (total_counts_tmp) ta <- add_total_counts(ta)
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$counts)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)
      
  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
     filter(taxon_id == spike_taxon) %>%
     select(sample_id, spike_count = count)
	        
  # calculate total density per sample
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(total_absolute_abundance = (spike_added * ( total_counts - spike_count ) / spike_count) / material_sampled)
			  
  # make counts table with total abundances
  tot_abs_abundance <- ta$samples %>%
     select(sample_id, total_absolute_abundance)
		        
  # calculate total for each taxon in each sample
  ta$counts <- ta$counts %>%
	  left_join(tot_abs_abundance)%>%
	  mutate(absolute_abundance = rel_abundance * total_absolute_abundance)
					  
  # remove spike_abundance
  ta$samples$spike_count <- NULL
					    
  # remove total_absolute_abundance from counts
  ta$counts$total_absolute_abundance <- NULL
				    
  # cleanup
  if (total_counts_tmp) ta$samples$total_counts <- NULL
  if (rel_abundance_tmp) ta$counts$rel_abundance <- NULL
					      
  # return ta object
  ta
					        
}
