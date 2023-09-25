#' Add sample table to the tidytacos object
#'
#' \code{add_sample_tibble} adds a sample tibble to the tidytacos object.
#'
#' This function adds a sample tibble containing metadata for each sample to the
#' tidytacos object. It is used after initiating a tidytacos object
#' using a numerical abundance matrix and the function
#' \code{\link{create_tidytacos}}. Also see \code{\link{add_taxon_tibble}}
#' to update the taxon data of the tidytacos object.
#'
#' @param ta tidytacos object.
#' @param sample_tibble A tibble containing sample data for each sample. samples
#'   should be rows, while sample data should be columns. At least one column
#'   name needs to be shared with the sample tibble of ta. The default shared
#'   column name is 'sample'.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Initiate sample tibble
#' sample <- c("sample1", "sample2")
#' environment <- c("food fermentation", "human stool")
#' sample_tibble <- tibble::tibble(sample, environment)
#'
#' # Add sample tibble to tidytacos object
#' data <- data %>%
#' add_sample_tibble(sample_tibble)
#'
#' @export
add_sample_tibble <- function(ta, sample_tibble) {

  purrr::modify_at(ta, "samples", left_join, sample_tibble)

}

#' Add total reads per sample
#'
#' \code{add_total_counts} adds the total reads per sample to the samples tibble
#' of a tidytacos object.
#'
#' This function adds the total reads per sample to the samples tibble of a
#' tidytacos object under the variable name total_counts.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total counts
#' data <- data %>%
#'  add_total_counts()
#'
#' @export
add_total_counts <- function(ta, step = "current") {

  # remove lib_size if already present
  ta$samples$total_counts <- NULL

  if (step == "current") {

    # make table with sample and library size
    lib_sizes <- ta$counts %>%
      group_by(sample_id) %>%
      summarize(total_counts = sum(count)) %>%
      select(sample_id, total_counts)

  } else {

    # make table with sample and library size
    step_of_interest <- step
    lib_sizes <- ta$lib_sizes %>%
      filter(step == step_of_interest) %>%
      select(sample_id, total_counts)

  }

  # add total counts to sample table
  ta$samples <-
    ta$samples %>%
    left_join(lib_sizes, by = "sample_id") %>%
    mutate(total_counts = ifelse(is.na(total_counts), 0, total_counts))

  # return ta object
  ta

}


#' Add alpha diversity measures
#'
#' \code{add_alphas} adds two alpha diversity measures to the
#' samples tibble of a tidytacos object.
#'
#' This function adds two alpha diversity measures (observed and inverse
#' Simpson) to the samples tibble of a tidytacos object under the variable
#' names observed and inverse_simpson, respectively. This function will also
#' add relative abundances if not present using \code{\link{add_rel_abundance}}.
#'
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_alphas()
#' @export
add_alphas <- function(ta) {

  # if rel abundances not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$counts)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # make table with sample, divObserved and divInvSimpson
  diversities <- ta$counts %>%
    filter(count > 0) %>%
    group_by(sample_id) %>%
    summarize(
      observed = n(),
      inverse_simpson = 1 / sum(rel_abundance ^ 2)
    ) %>%
    ungroup()

  # add diversity measure to sample table
  ta$samples = left_join(ta$samples, diversities, by = "sample_id")

  # cleanup
  if (rel_abundance_tmp) ta$counts$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add clustered sample order
#'
#' \code{add_sample_clustered} adds a new variable defining a sample order based
#' on similarity after clustering to the samples tibble of a tidytacos
#' object.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' hierarchical average linkage clustering of samples. It will then add a new
#' factor variable "samples_clustered" to the samples tibble of a tidytacos
#' object. This function is extremely useful if one wants to plot similar
#' samples together.
#'
#' @importFrom stats hclust
#' @param ta tidytacos object.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_sample_clustered()
#'
#' @export
add_sample_clustered <- function(ta) {

  # make relative abundance matrix
  rel_abundance_matrix <- rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix <- vegdist(rel_abundance_matrix, method = "bray")

  # perform hierarchical clustering
  clust <- hclust(dist_matrix, method = "average")

  # make table with samples in order of clustering
  samples_clustered <- tibble(
    sample_id = clust$labels[clust$order],
    sample_clustered = factor(sample_id, levels = sample_id)
  )

  # add sample_clustered to samples table
  ta$samples <- ta$samples %>%
    left_join(samples_clustered, by = "sample_id")

  # return ta object
  ta

}

# Helper function to prepare 2/3D coordinates of ord
get_dimensions <- function(dim_df, names, dims) {

  ordnames <- c("ord1", "ord2")
  if (dims == 3) {
    ordnames <- c(ordnames, "ord3")
  }

  dim_df %>%
    `colnames<-`(ordnames) %>%
    as_tibble() %>%
    mutate(sample_id = !! names)
}
# Calculate pcoa coordinates and variances
perform_pcoa <- function(ta, dist_matrix, dims=2, ...){

  ord <- list()
  pcoa <- stats::cmdscale(dist_matrix, k = dims, eig = T, list = T, ...)
  ord$variances <- pcoa$eig / sum(pcoa$eig)
  ord$dimensions <- get_dimensions(
    pcoa$points, rownames(pcoa$points), dims=dims)
  ord
}

# Calculate tsne coordinates and variances
perform_tsne <- function(ta, dist_matrix, dims=2, ...) {
  force_optional_dependency("Rtsne")

  ord <- list()
  tsne <- Rtsne::Rtsne(dist_matrix, dims=dims, ...)
  ord$dimensions <- get_dimensions(
    tsne$Y, rownames(as.matrix(dist_matrix)), dims = dims)
  ord$variances <- tsne$costs / sum(tsne$costs)
  ord
}

# Calculate umap coordinates and variances
perform_umap <- function(ta, dist_matrix, dims=2, ...) {
  force_optional_dependency("umap")
  ord <- list()
  umap <- umap::umap(as.matrix(dist_matrix), n_components=dims, ...)
  ord$dimensions <- get_dimensions(
    umap$layout, rownames(umap$layout), dims = dims)
  ord$variances <- umap$knn$distances / sum(umap$knn$distances)
  ord
}

#' Add dimensionality ordination
#'
#' \code{add_ord} adds the first x dimensions of a dimensionality reduction method
#' of a given dissimilarity matrix to two new variables of the
#' samples tibble of a tidytacos object.
#'
#' This function calculates the distance between samples followed by
#' an ordination analysis. It will then add the first n dimensions to
#' the samples tibble of a tidytacos object named "ord1", "ord2", ... This
#' function will also add relative abundances if not present using
#' \code{\link{add_rel_abundance}}.
#' @param ta tidytacos object.
#' @param distance the distance indices to use, see \code{\link[vegan]{vegdist}} 
#' @param method the ordination method to use to calculate coordinates. Choice from pcoa, tsne, umap
#' @param dims the amount of dimensions to reduce the distances to.
#' @param binary perform presence/absence standardisation before distance computation.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356, 456, 678),
#'  ncol = 3
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2", "sample3")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add pcoa
#' data <- data %>%
#'  add_ord()
#'
#' @export
add_ord <- function(ta, distance="bray", method="pcoa", dims=2, binary=FALSE, ...) {

  methods = c("pcoa", "tsne", "umap")
  if (!method %in% methods) {
    stop(paste("Select a method from", paste0(method, collapse=",")))
  }

  # if add_ord was run before, remove coordinates from sample table
  if ("ord_method" %in% names(ta)) {
    warning("Overwriting previous ord data")
    ta$samples <- ta$samples %>% 
        select(-num_range("ord", 0:length(ta$samples$sample_id)))
  }

  # make relative abundance matrix
  rel_abundance_matrix <- rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix = vegan::vegdist(rel_abundance_matrix, method = distance, binary=binary)

  if (method == "pcoa") {
    ord <- perform_pcoa(ta, dist_matrix, dims=dims, ...)
  }

  if (method == "tsne") {
    ord <- perform_tsne(ta, dist_matrix, dims=dims, ...)
  }

  if (method == "umap") {
    ord <- perform_umap(ta, dist_matrix, dims=dims, ...)
  }

  # add ord dimensions to sample table
  ta$samples <- ta$samples %>%
    left_join(ord$dimensions, by = "sample_id")

  # add ord variances to ta object
  ta$ord_variances <- ord$variances
  ta$ord_method <- method
  # return ta object
  ta

}

#' Add spike ratio
#'
#' \code{add_spike_ratio} adds a new variable showing the ratio total counts
#' to spike counts to the samples tibble of a tidytacos object.
#'
#' This function calculates the spike ratio defined as the total sample
#' counts to the spike counts and adds this as a new variable
#' "spike_ratio" to the samples tibble of a tidytacos object. This function
#' is useful if a DNA spike was added prior to sequencing and is based on the
#' method described by
#' \href{https://doi.org/10.1016/j.soilbio.2016.02.003}{Smets et al., 2016}.
#'
#' Credits to Wenke Smets for the idea of spiking samples prior to 16S
#' sequencing and the initial implementation of this function.
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon_id of the spike.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_spike_ratio(spike_taxon = "t1")
#
#' @export
add_spike_ratio <- function(ta, spike_taxon) {

  # if lib_size not present: add temporarily
  lib_size_tmp <- ! "total_counts" %in% names(ta$samples)
  if (lib_size_tmp) ta <- add_total_counts(ta)

  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_abundance = count)

  # calculate spike ratio (non-spike abundance to spike abundance)
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(spike_ratio = ( total_counts - spike_abundance ) / spike_abundance)

  # remove spike_abundance
  ta$samples$spike_abundance <- NULL

  # cleanup
  if (lib_size_tmp) ta$samples$total_counts <- NULL

  # return ta object
  ta

}

#' Add cluster number
#'
#' \code{add_sample_cluster} adds a new variable to the samples tibble of a
#' tidytacos object defining to what cluster a sample belongs.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' hierarchical average linkage clustering of samples. The user provides a
#' number of desired clusters which will be used to assign the samples to. A new
#' variable named "cluster" will be added to the samples tibble of a
#' tidytacos object defining to what cluster a sample belongs.
#'
#' @param ta tidytacos object.
#' @param n_clusters Numerical. Number of desired clusters.
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
#'                      taxa_are_columns = FALSE
#'                      )
#'
#' # Add total abundance
#' data <- data %>%
#'  add_sample_cluster(n_clusters = 2)
#'
# Adds a variable "cluster" to the samples table
# To do: merge with add_sample_clustered somehow
#
#' @importFrom stats cutree
#' @export
add_sample_cluster <- function(ta, n_clusters) {

  # make relative abundance matrix
  rel_abundance_matrix <- rel_abundance_matrix(ta)

  # make Bray-Curtis distance matrix
  dist_matrix <- vegdist(rel_abundance_matrix, method = "bray")

  # perform hierarchical clustering
  clust <- hclust(dist_matrix, method = "average")

  samples_clusters <-
    tibble(
      sample_id = clust$labels,
      cluster = cutree(clust, k = n_clusters)
    ) %>%
    mutate(cluster = str_c("cluster", cluster, sep = " "))

  ta$samples <-
    left_join(ta$samples, samples_clusters, by = "sample_id")

  ta

}


#' Add total absolute abundances to samples table
#'
#' \code{add_total_absolute_abundance} adds total absolute abundance to the samples table of a
#' tidytacos object.
#'
#' This function adds the total absolute abundances to the samples table
#' of a tidytacos object under the variable name "total_absolute_abundance".
#'
#' @param ta tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the samples table which indicates how much spike was added per sample, e.g. 16S rRNA gene copy numbers added to the DNA extraction tube.
#'
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'   c(1500, 1300, 14, 280, 356, 9),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2", "taxon3")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'   taxa_are_columns = FALSE
#' )
#' data$samples$spike_added <- c(100, 150)
#'
#' # Add total abundance
#' data <- data %>%
#'   add_total_absolute_abundance(spike_taxon = "t3")
#'
#' @export
add_total_absolute_abundance <- function(ta, spike_taxon, spike_added = spike_added) {
  spike_added <- rlang::enquo(spike_added)
  
  if (!rlang::quo_name(spike_added) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::quo_name(spike_added),
      "that defines the quantity of spike added to the sample."
    ))
  }
  
  # if total_counts not present: add temporarily
  total_counts_tmp <- !"total_counts" %in% names(ta$samples)
  if (total_counts_tmp) ta <- add_total_counts(ta)
  
  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_count = count)
  
  # calculate total absolute abundance per sample
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(total_absolute_abundance = (!!spike_added * (total_counts - spike_count) / spike_count))
  
  # remove spike_abundance
  ta$samples$spike_count <- NULL
  
  # cleanup
  if (total_counts_tmp) ta$samples$total_counts <- NULL
  
  # Warn about samples without spike
  samples_w_no_spike <- unique(ta$samples$sample_id[which(is.na(ta$samples$total_absolute_abundance))])
  if (length(samples_w_no_spike) > 0) {
    warning(
      paste(
        "Sample without spike taxon detected:",
        format(samples_w_no_spike, trim = TRUE), "\n"
      )
    )
  }
  
  # return ta object
  ta
}


#' Add total densities to samples table
#'
#' \code{add_total_density} adds total density to the samples table of a
#' tidytacos object.
#'
#' This function adds the total densities to the samples table
#' of a tidytacos object under the variable name "total_density".
#'
#' @param ta tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the samples table which indicates how much spike was added per sample, e.g. 16S rRNA gene copy numbers added to the DNA extraction tube.
#' @param material_sampled The column name indicating the amount of material from which DNA was extracted, e.g gram of soil. This parameter encourages researchers to consider that absolute abundances are only meaningful if they can be translated into densities.
#'
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'   c(1500, 1300, 14, 280, 356, 9),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2", "taxon3")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'   taxa_are_columns = FALSE
#' )
#' data$samples$spike_added <- c(100, 150)
#' data$samples$material_sampled <- c(1, 5)
#'
#' # Add total abundance
#' data <- data %>%
#'   add_total_density(spike_taxon = "t3")
#'
#' @export
add_total_density <- function(ta, spike_taxon, spike_added = spike_added, material_sampled = material_sampled) {
  spike_added <- rlang::enquo(spike_added)
  material_sampled <- rlang::enquo(material_sampled)
  
  if (!rlang::quo_name(spike_added) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::quo_name(spike_added),
      "that defines the quantity of spike added to the sample."
    ))
  }
  
  if (!rlang::quo_name(material_sampled) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::quo_name(material_sampled),
      "that defines the quantity of sample used."
    ))
  }
  
  # if total_counts not present: add temporarily
  total_counts_tmp <- !"total_counts" %in% names(ta$samples)
  if (total_counts_tmp) ta <- add_total_counts(ta)
  
  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_count = count)
  
  # calculate total absolute abundance per sample
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(total_density = (!!spike_added * (total_counts - spike_count) / spike_count)/ !!material_sampled)
  
  # remove spike_abundance
  ta$samples$spike_count <- NULL
  
  # cleanup
  if (total_counts_tmp) ta$samples$total_counts <- NULL
  
  # Warn about samples without spike
  samples_w_no_spike <- unique(ta$samples$sample_id[which(is.na(ta$samples$total_density))])
  if (length(samples_w_no_spike) > 0) {
    warning(
      paste(
        "Sample without spike taxon detected:",
        format(samples_w_no_spike, trim = TRUE), "\n"
      )
    )
  }
  
  # return ta object
  ta
}


