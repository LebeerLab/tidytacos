#' Add metadata to the tidytacos object
#'
#' \code{add_metadata} adds sample or taxon metadata to the sample or taxon
#' table, respectively, of a tidytacos object.
#'
#' @param ta A tidytacos object.
#' @param metadata A tibble containing data for each sample or taxon.
#'   Samples/taxa should be rows, while metadata variables should be columns. At
#'   least one column name needs to be shared with the sample or taxa table of
#'   the tidytacos object. The default shared column name is 'sample' for
#'   samples and 'taxon' for taxa.
#' @param table_type The type of table to add, either 'sample' or 'taxa'.
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
#' add_metadata(sample_tibble)
#'
#' @export
add_metadata <- function(ta, metadata_tibble, table_type = "sample") {

  if (table_type == "sample") {
    purrr::modify_at(ta, "samples", left_join, metadata_tibble)
  } else if (table_type == "taxa") {
    ta <- purrr::modify_at(ta, "taxa", left_join, metadata_tibble)
    infer_rank_names(ta)
  } else {
    stop("table_type must be either 'sample' or 'taxa'")
  }

}

#' Add total read count per sample
#'
#' \code{add_total_count} adds the total read count per sample to the sample
#' table of a tidytacos object under the variable name total_count.
#'
#' @param ta A tidytacos object.
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
#'  add_total_count()
#'
#' @export
add_total_count <- function(ta) {

  # make table with sample and total count
  lib_sizes <- ta$counts %>%
    group_by(sample_id) %>%
    summarize(total_count = sum(count)) %>%
    select(sample_id, total_count)

  # add total count to sample table
  ta$samples <-
    ta$samples %>%
    left_join(lib_sizes, by = "sample_id") %>%
    mutate(total_count = ifelse(is.na(total_count), 0, total_count))

  # return ta object
  ta

}

#' Add alpha diversity measure
#'
#' \code{add_alpha} adds an alpha diversity measures to the sample table of a
#' tidytacos object.
#'
#' This function can add different alpha diversity measures to the sample table, specified by the method argument.
#' The following methods are available:
#' - invsimpson: Inverse Simpson index
#' - shannon: Shannon index
#' - simpson: Simpson index
#' - pielou: Pielou's evenness index
#' - obs: Observed richness
#' - s.chao1: Chao1 richness estimator
#' - s.ace: ACE richness estimator
#'
#' @param ta A tidytacos object.
#' @param method The diversity measure to use, see \code{\link[vegan]{diversity}} for further information on these.
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
#'  add_alpha()
#'
#' @export
add_alpha <- function(ta, method="invsimpson") {

  vegan_standard_methods <- c("invsimpson","shannon","simpson")
  vegan_estimateR_methods <- c("obs", "s.chao1", "s.ace")

  method <- tolower(method)

  if (!method %in% lapply(alpha_metrics, tolower)) {
    stop(paste("Select a method from", paste0(alpha_metrics, collapse=", ")))
  }

  if (method %in% vegan_standard_methods) {
    M <- ta %>% counts_matrix(sample_name=sample_id, taxon_name=taxon_id)
    D <- vegan::diversity(M, index=method)
    diversities <- tibble::tibble(sample_id = names(D), !!method := D)
  }

  if (method %in% vegan_estimateR_methods) {
    M <- ta %>% counts_matrix(sample_name=sample_id, taxon_name=taxon_id)
    D <- vegan::estimateR(M)
    selection <- grepl(method, rownames(D), ignore.case = TRUE)
    selection_names <- rownames(D)[selection]
    Dt <- D %>% t()
    Dt <- Dt[,selection_names]

    diversities <- Dt %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(sample_id = names(D[1,]))

    if (is.null(dim(Dt))) {
      diversities <- diversities %>% 
        dplyr::rename(!!method := value)
    }
  }

  if (method == "pielou") diversities <- calculate_alpha_pielou(ta)

  # add diversity measure to sample table
  ta$samples = left_join(ta$samples, diversities, by = "sample_id")

  # return ta object
  ta

}

#' Add alpha diversity measures
#'
#' \code{add_alpha} adds selected alpha diversity measures to the sample table of a
#' tidytacos object.
#'
#' This function can add multiple different alpha diversity measures to the sample table, specified by the methods argument.
#' @param ta A tidytacos object.
#' @param methods A character vector of the diversity measure to use, see \code{\link[tidytacos]{add_alpha}} for examples.
#' Optionally use 'all' to add all diversity measures.
#' @export
add_alphas <- function(ta, methods="all") {

  if (any(methods == "all")) {
    methods <- alpha_metrics
  }

  for (method in methods) {
    ta <- add_alpha(ta, method)
  }
  ta
}


calculate_alpha_pielou <- function(ta) {

  M <- ta %>% counts_matrix(sample_name=sample_id, taxon_name=taxon_id)
  H <- vegan::diversity(M, index="shannon")
  J <- H / log(vegan::specnumber(M))

  tibble(sample_id = names(J), pielou = J)
}

#' Add clustering-based sample order
#'
#' \code{add_sample_clustered} adds a new variable defining a sample order based
#' on a hierarchical clustering of the samples.
#'
#' This function calculates the Bray-Curtis distances between samples followed
#' by hierarchical average linkage clustering of samples. It will then add a new
#' factor variable "sample_clustered" to the sample tibble of a tidytacos
#' object. This function is useful if one wants to plot similar samples
#' together.
#'
#' @param ta  A tidytacos object.
#'
#' @importFrom stats hclust
#' @export
add_sample_clustered <- function(ta) {

  # if only one sample => no clustering
  if (length(ta$samples$sample_id) == 1) {
    ta$samples$sample_clustered <- factor(ta$samples$sample_id)
    return(ta) 
  } 

  # make relative abundance matrix
  rel_abundance_matrix <- rel_abundance_matrix(ta, sample_name=sample_id, taxon_name=taxon_id)

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

# Helper function to prepare 2/3D coordinates of ord.
get_dimensions <- function(dim_df, names, dims) {

  ordnames <- c("ord1", "ord2")

  if (dims >= 3) {
    for (i in 3:dims){
      ordnames <- c(ordnames, paste0("ord", i))
    }
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

#' Add ordination
#'
#' \code{add_ord} adds the first n dimensions of a dimensionality reduction
#' method performed on a given dissimilarity matrix as new variables to the
#' sample table of a tidytacos object.
#'
#' This function calculates the dissimilarities between samples followed by
#' an ordination analysis. It will then add the first n dimensions to
#' the sample table of a tidytacos object named "ord1", "ord2", ... This
#' function will also add relative abundances if not present using
#' \code{\link{add_rel_abundance}}.
#'
#' @param ta A tidytacos object.
#' @param distance The distance indices to use, see
#'   \code{\link[vegan]{vegdist}}.
#' @param method The ordination method to use to calculate coordinates. Choices
#'   are "pcoa", "tsne", "umap".
#' @param dims The amount of dimensions to reduce the dissimilarities to.
#' @param binary Perform presence/absence standardisation before distance
#'   computation or not.
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
  rel_abundance_matrix <- rel_abundance_matrix(ta, sample_name=sample_id, taxon_name=taxon_id)

  if (distance == "aitchison") {
    # Euclidean distance between CLR transformed abundances
    rel_abundance_matrix <- rel_abundance_matrix %>%
      vegan::decostand(method = "clr", pseudocount = 1)
    dist_matrix <- vegan::vegdist(rel_abundance_matrix, method = "euclidean")
  } else {
    # make distance matrix
    dist_matrix = vegan::vegdist(rel_abundance_matrix, method = distance, binary=binary)
  }

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
#' \code{add_spike_ratio} calculates the ratio of non-spike to spike reads for
#' each sample and adds this to the sample table under the name "spike_ratio".
#'
#' This function is useful if a DNA spike was added prior to sequencing and is
#' based on the method described by
#' \href{https://doi.org/10.1016/j.soilbio.2016.02.003}{Smets et al., 2016}.
#'
#' Without calculating absolute abundances, the spike ratio allows to compare absolute abundances between sample. For example, if the spike ration of one sample is twice that of another, then the absolute number of sequenced strands at the time of spiking in the one sample is twice that of the other sample.
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
  lib_size_tmp <- ! "total_count" %in% names(ta$samples)
  if (lib_size_tmp) ta <- add_total_count(ta)

  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_abundance = count)

  # calculate spike ratio (non-spike abundance to spike abundance)
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(spike_ratio = (total_count - spike_abundance) / spike_abundance)

  # remove spike_abundance
  ta$samples$spike_abundance <- NULL

  # cleanup
  if (lib_size_tmp) ta$samples$total_count <- NULL

  # return ta object
  ta

}

#' Clusters samples into n clusters
#'
#' \code{cluster_samples} clusters the samples into n clusters and adds these
#' clusters to a new variable "cluster" in the sample table.
#'
#' This function calculates the Bray-Curtis distance between samples followed by
#' hierarchical average linkage clustering of samples. The user provides a
#' number of desired clusters which will be used to assign the samples to. A new
#' variable named "cluster" will be added to the samples tibble of a
#' tidytacos object defining to what cluster a sample belongs.
#'
#' @param ta A tidytacos object.
#' @param n_clusters The number of desired clusters.
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
#'  cluster_samples(n_clusters = 2)
#'
# Adds a variable "cluster" to the sample table
# To do: merge with add_sample_clustered somehow
#
#' @importFrom stats cutree
#' @export
cluster_samples<- function(ta, n_clusters) {

  # make relative abundance matrix
  rel_abundance_matrix <- rel_abundance_matrix(ta, sample_name=sample_id, taxon_name=taxon_id)

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


#' Add total absolute abundances of samples
#'
#' \code{add_total_absolute_abundance} calculates the total absolute abundances
#' of the samples given a spike taxon, and adds this to the sample table under
#' the column name "total_absolute_abundance".
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the samples table which indicates how
#'   much spike was added per sample, e.g. 16S rRNA gene copy numbers added to
#'   the DNA extraction tube.
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

  # if total_count not present: add temporarily
  total_count_tmp <- !"total_count" %in% names(ta$samples)
  if (total_count_tmp) ta <- add_total_count(ta)

  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_count = count)

  # calculate total absolute abundance per sample
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(total_absolute_abundance = (!!spike_added * (total_count - spike_count) / spike_count))

  # remove spike_abundance
  ta$samples$spike_count <- NULL

  # cleanup
  if (total_count_tmp) ta$samples$total_count <- NULL

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


#' Add total densities of samples
#'
#' \code{add_total_density} adds the total microbial density to the sample table
#' of a tidytacos object under the column name "total_density".
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the samples table which indicates how
#'   much spike was added per sample, e.g. 16S rRNA gene copy numbers added to
#'   the DNA extraction tube.
#' @param material_sampled The column name indicating the amount of material
#'   from which DNA was extracted, e.g gram of soil. This parameter encourages
#'   researchers to consider that absolute abundances are only meaningful if
#'   they can be translated into densities.
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

  # if total_count not present: add temporarily
  total_count_tmp <- !"total_count" %in% names(ta$samples)
  if (total_count_tmp) ta <- add_total_count(ta)

  # make sample table with spike abundances
  spike_counts <- ta$counts %>%
    filter(taxon_id == spike_taxon) %>%
    select(sample_id, spike_count = count)

  # calculate total absolute abundance per sample
  ta$samples <- ta$samples %>%
    left_join(spike_counts, by = "sample_id") %>%
    mutate(total_density = (!!spike_added * (total_count - spike_count) / spike_count)/ !!material_sampled)

  # remove spike_abundance
  ta$samples$spike_count <- NULL

  # cleanup
  if (total_count_tmp) ta$samples$total_count <- NULL

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

#' Perform anosim test
#'
#' \code{perform_anosim} performs the anosim test for statistical difference
#' between groups of samples. The null hypothesis is that there is no difference
#' between microbial communities in the groups of samples.
#'
#' @param ta A tidytacos object.
#' @param group A column in the sample table to group the samples on.
#' @param distance The dissimilarity measure to use.
#' @param permutations The number of permutations.
#'
#' @export
perform_anosim <- function(ta, group, ...){

  M <- ta %>% counts_matrix()
  group <- rlang::enquo(group)

  if (length(M[,1]) < length(ta$samples$sample_id)) {
        warning("Empty samples found, ignoring them in analysis")
        ta <- ta %>% remove_empty_samples()
  }

  vegan::anosim(M, ta$samples %>% pull(!!group), ...)

}
