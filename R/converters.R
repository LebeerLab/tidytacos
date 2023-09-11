#' Initiate tidytacos object
#'
#' \code{tidytacos} returns a tidytacos object given a numeric matrix.
#'
#' This function initiates a tidytacos object based on a numeric matrix. It
#' will automatically create a dummy taxa and sample table which will need to be
#' updated using the functions \code{\link{add_taxon_tibble}} and
#' \code{\link{add_sample_tibble}}.
#'
#' @param counts_matrix Numerical matrix containing the count data.
#' @param taxa_are_columns A logical scalar. Are the taxa defined in columns?
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
#'
#' \dontrun{
#' tidytacos("a")
#' }
#'
#' @export
create_tidytacos <- function(counts_matrix, taxa_are_columns = TRUE) {

  if (
    ! is.matrix(counts_matrix) |
    ! is.numeric(counts_matrix)
  ) stop("first argument should be a numeric matrix")

  if (! taxa_are_columns) counts_matrix = t(counts_matrix)

  counts_matrix <-
    counts_matrix[, colSums(counts_matrix) != 0]

  ta <- list()
  class(ta) <- "tidytacos"

  n_samples <- nrow(counts_matrix)
  n_taxa <- ncol(counts_matrix)

  sample_ids <- str_c("s", 1:n_samples)
  taxon_ids <- str_c("t", 1:n_taxa)

  ta$counts <-
    counts_matrix %>%
    as.vector() %>%
    tibble(count = .) %>%
    mutate(sample_id = rep(!! sample_ids, times = !! n_taxa)) %>%
    mutate(taxon_id = rep(!! taxon_ids, each = !! n_samples)) %>%
    filter(count > 0)

  ta$samples <-
    counts_matrix %>%
    rownames() %>%
    tibble(sample = .) %>%
    mutate(sample_id = !! sample_ids)

  ta$taxa <-
    counts_matrix %>%
    colnames() %>%
    tibble(taxon = .) %>%
    mutate(taxon_id = !! taxon_ids)

  ta

}

#' Write community data in tidytacos format
#' @importFrom readr write_csv
#' @param ta a tidytacos object
#' @param dout the directory to store the three tidytacos tables in
#' @export
write_tidytacos <- function(ta, dout) {
  if (!dir.exists(dout)) {dir.create(dout)}
  write_csv(ta$samples, paste0(dout, "/samples.csv"))
  write_csv(ta$taxa, paste0(dout, "/taxa.csv"))
  write_csv(ta$counts, paste0(dout, "/counts.csv"))
}

#' Read community data written by tidytacos
#' @importFrom readr read_csv
#' @param din directory containing the a sample, taxa and counts table in csv format
#' @param samples the name of the samples table, defaults to samples.csv
#' @param taxa the name of the taxa table, defaults to taxa.csv
#' @param counts the name of the counts table, defaults to counts.csv
#' @export
read_tidytacos <- function(din, samples = "samples.csv", taxa = "taxa.csv",
                               counts = "counts.csv") {
  samples <- readr::read_csv(paste0(din, "/", samples), col_types = readr::cols())
  taxa <- readr::read_csv(paste0(din, "/", taxa), col_types = readr::cols())
  counts <- readr::read_csv(paste0(din, "/", counts), col_types = readr::cols())
  make_tidytacos(
    samples, taxa, counts, sample_name = sample_id, taxon_name = taxon_id
  )
}

#' Reset the taxon and sample IDs
#' @param ta a tidytacos object
#' @export
reset_ids <- function(ta, keep_prev = F) {

  if (keep_prev) {

    ta <-
      ta %>%
      mutate_samples(sample_id_prev = sample_id) %>%
      mutate_taxa(taxon_id_prev = taxon_id)

  }

  ta %>%
    change_id_samples(sample_id_new = str_c("s", seq_len(n()))) %>%
    change_id_taxa(taxon_id_new = str_c("t", seq_len(n())))

}

#' Convert tidytacos object to phyloseq object
#'
#' \code{as_phyloseq} returns a phyloseq object given a tidytacos object.
#'
#' This function will convert a tidytacos object into a phyloseq object for
#' alternative processing using the phyloseq package. To convert from a phyloseq
#' object to a tidytacos object use \code{\link{as_tidytacos}}.
#'
#' @param ta tidytacos object.
#' @param sample  The sample names required for a phyloseq object. Default is
#'   "sample" column in sample tibble of the tidytacos object.
#' @param taxon The taxon names required for a phyloseq object. Default is
#'   "taxon" column in taxon tibble of the tidytacos object.
#'
#' @export
as_phyloseq <- function(ta, sample = sample, taxon = taxon) {

  force_optional_dependency("phyloseq")
  if ("phyloseq" %in% class(ta)) return(ta)

  sample <- rlang::enexpr(sample)
  taxon <- rlang::enexpr(taxon)

  ta <- change_id_samples(ta, sample_id_new = !! sample)
  ta <- change_id_taxa(ta, taxon_id_new = !! taxon)

  otu_table <-
    ta$counts %>%
    spread(key = taxon_id, value = count, fill = 0) %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$sample_id) %>%
    select(- sample_id) %>%
    as.matrix() %>%
    phyloseq::otu_table(taxa_are_rows = F)

  if (ncol(ta$samples) == 1) {
    ta <- mutate_samples(ta, dummy = as.character(1:nrow(ta$samples)))
  }

  if (ncol(ta$taxa) == 1) {
    ta <- mutate_taxa(ta, dummy = as.character(1:nrow(ta$taxa)))
  }

  sample_data <-
    ta$samples %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$sample_id) %>%
    select(- sample_id) %>%
    phyloseq::sample_data()

  tax_table <-
    ta$taxa %>%
    `attr<-`("class", "data.frame") %>%
    `rownames<-`(.$taxon_id) %>%
    select(- taxon_id) %>%
    as.matrix() %>%
    phyloseq::tax_table()

  phyloseq::phyloseq(otu_table, sample_data, tax_table)

}

#' Convert phyloseq object to tidytacos object
#'
#' \code{as_tidytacos} returns a tidytacos object given a phyloseq
#' object.
#'
#' This function will convert a phyloseq object into a tidytacos object. To
#' convert from a tidytacos object to a phyloseq object use
#' \code{\link{as_phyloseq}}.
#'
#' @param ps Phyloseq object.
#'
#' @export
as_tidytacos <- function(ps) {

  if ("tidytacos" %in% class(ps)) return(ps)

  # convert sample data to tibble
  samples <-
    phyloseq::sample_data(ps)@.Data %>%
    `names<-`(phyloseq::sample_data(ps)@names) %>%
    do.call(what = tibble) %>%
    mutate(sample = phyloseq::sample_data(ps)@row.names)

  # convert taxon table to tibble
  taxa <-
    phyloseq::tax_table(ps)@.Data %>%
    as_tibble() %>%
    mutate(taxon = phyloseq::tax_table(ps) %>% row.names()) %>%
    `names<-`(names(.) %>% str_to_lower())

  # make sure that taxa are columns in counts table
  if (phyloseq::taxa_are_rows(ps)) {
    phyloseq::otu_table(ps) <- phyloseq::t(phyloseq::otu_table(ps))
  }

  phyloseq::otu_table(ps)@.Data %>%
    create_tidytacos() %>%
    add_sample_tibble(samples) %>%
    add_taxon_tibble(taxa)

}

#' Convert matrix with counts to tidy data frame
#'
#' \code{as_counts} returns a tidy data frame given a numerical counts
#' matrix.
#'
#' This function will convert a numerical counts matrix into a tidy data
#' frame. To convert a tidy data frame into a numerical counts matrix
#' use \code{\link{as_counts_matrix}}.
#'
#' @param counts_matrix The ambundance matrix that will be converted.
#' @param taxa_are_columns A logical scalar. Are the taxa defined in columns?
#'   Default is TRUE.
#' @param value Name of resulting colum containing the count data. Default
#'   is "counts".
#'
#' @export
as_counts <- function(counts_matrix, taxa_are_columns = TRUE,
                          value = "counts") {

  if (
    ! is.matrix(counts_matrix) |
    ! is.numeric(counts_matrix)
  ) stop("first argument should be an counts matrix")

  if (! taxa_are_columns) counts_matrix = t(counts_matrix)

  counts_matrix %>%
    as_tibble() %>%
    mutate(sample_id = row.names(counts_matrix)) %>%
    gather(key = "taxon_id", value = !! value, - sample_id) %>%
    filter(!! value > 0)

}

#' Convert counts tidy data frame to matrix.
#'
#' \code{as_counts_matrix} returns a numerical matrix given a tidy
#' counts data frame.
#'
#' This function will convert a counts tidy data frame into a numerlical
#' counts matrix. To convert a numerical counts matrix into a counts
#' tidy data frame use \code{\link{as_counts_matrix}}.
#'
#' @param counts The counts tidy data frame that will be converted.
#' @param value Name of colum containing the counts data. Default is
#'   "counts".
#'
#' @export
as_counts_matrix <- function(counts, value = count) {

  if (
    ! is.data.frame(counts) |
    is.null(counts$taxon_id) |
    is.null(counts$sample_id)
  ) stop("first argument should be a counts table (data frame)")

  value <- enquo(value)

  counts_wide <- counts %>%
    select(sample_id, taxon_id, !! value) %>%
    spread(key = taxon_id, value = !! value, fill = 0)

  counts_wide %>%
    select(- sample_id) %>%
    as.matrix() %>%
    `row.names<-`(counts_wide$sample_id)

}

#' Merge two tidytacos objects.
#'
#' \code{merge_tidytacos} merges two tidytacos objects and returns one
#' single tidytacos object.
#'
#' This function will merge two tidytacos objects into one. It is useful if
#' one wants to merge data obtained from different sequencing runs. Therefore,
#' this function requirers that both tidytacos objects contain a "run"
#' variable in their samples table, indicating their origin.
#'
#' @param ta1 The first tidytacos object.
#' @param ta2 The second tidytacos object.
#'
#' @export
merge_tidytacos <- function(ta1, ta2, taxon_identifier = sequence) {

  taxon_identifier <- rlang::ensym(taxon_identifier)

  ti <- rlang::as_string(taxon_identifier)
  if (! (ti %in% names(ta1$taxa) & ti %in% names(ta2$taxa))) {
    stop("the taxon identifier was not found in one or both of the ta objects")
  }

  # make sure that sample names are unique
  ta1 <- change_id_samples(ta1, paste("ta1", sample_id, sep = "_"))
  ta2 <- change_id_samples(ta2, paste("ta2", sample_id, sep = "_"))

  # merge sample tables
  samples <- bind_rows(ta1$samples, ta2$samples)

  # change taxon ids to something meaningful across ta objects
  ta1 <- change_id_taxa(ta1, taxon_id_new = !! taxon_identifier)
  ta2 <- change_id_taxa(ta2, taxon_id_new = !! taxon_identifier)

  # merge taxa tables
  taxa <-
    bind_rows(ta1$taxa, ta2$taxa) %>%
    group_by(taxon_id) %>%
    summarize_all(function(x) {
      x <- unique(x)
      x <- x[! is.na(x)]
      if (length(x) == 1) return(x)
      as.character(NA)
    })

  # merge counts tables
  counts <- bind_rows(ta1$counts, ta2$counts)

  # make new ta object
  ta <- list(samples = samples, taxa = taxa, counts = counts)
  class(ta) <- "tidytacos"

  # give new sample names in new ta object
  ta <- reset_ids(ta)

  # return ta object
  ta

}

#' Create a tidytacos object from three tidy tables
#'
#' @param samples A tidy table containing sample information
#' @param taxa A tidy table containing taxon information
#' @param counts A tidy table, where each row represents the counts of a taxon in a sample
#' @param sample_name The column in the sample table that contains a unique identifier for each sample
#' @param taxon_name The column in the taxon table that contains a unique identifier for each taxon
#' @export
make_tidytacos <- function(samples, taxa, counts,
                               sample_name = sample, taxon_name = taxon) {

  sample_name <- rlang::enexpr(sample_name)
  taxon_name <- rlang::enexpr(taxon_name)

  list(samples = samples, taxa = taxa, counts = counts) %>%
    purrr::modify_at("samples", mutate, sample_id = !! sample_name) %>%
    purrr::modify_at("taxa", mutate, taxon_id = !! taxon_name) %>%
    purrr::modify_at("counts", rename, sample_id = !! sample_name) %>%
    purrr::modify_at("counts", rename, taxon_id = !! taxon_name) %>%
    purrr::modify_at("counts", filter, count > 0) %>%
    purrr::modify_at("counts", filter, sample_id %in% .$samples$sample_id) %>%
    purrr::modify_at("counts", filter, taxon_id %in% .$taxa$taxon_id) %>%
    `class<-`("tidytacos") %>%
    reset_ids()

}

