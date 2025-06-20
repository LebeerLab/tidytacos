#' Return some descriptive numbers
#'
#' `tacosum()` returns the number of samples, taxa and reads in the
#' tidytacos object.
#'
#' @param ta A tidytacos object.
#'
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'   c(1500, 1300, 280, 356),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x, taxa_are_columns = FALSE)
#'
#' # Report numbers
#' numbers <- data %>% tacosum()
#'
#' @export
tacosum <- function(ta) {
  c(
    n_samples = nrow(ta$samples),
    n_taxa = nrow(ta$taxa),
    n_reads = sum(ta$counts$count)
  )
}

#' Get beta diversity table
#'
#' `betas()` returns a tidy table with the beta diversity for each
#' combination of samples.
#'
#' This function calculates the beta diversity using the
#' [vegan::vegdist()] function of Vegan. It will report one diversity
#' estimate for each combination of samples.
#'
#'
#' @param ta A tidytacos object.
#' @param unique A logical scalar. Avoid redundancy by removing all self sample
#'   comparisons and keep only one of two pairwise comparisons? Default is TRUE.
#' @param method The dissimilarity index. See [vegan::vegdist()] for
#'   all options. Default is "bray".
#' @param binary A logical scalar. Perform presence/absence standardization
#'   before analysis. See [vegan::vegdist()]. Default is FALSE.
#' @inheritDotParams vegan::vegdist
#' @examples
#' # Initiate counts matrix
#' x <- matrix(
#'   c(1500, 1300, 280, 356),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <-
#'   create_tidytacos(x, taxa_are_columns = FALSE)
#'
#' # Report numbers
#' numbers <- data %>% betas()
#'
#' @export
betas <- function(ta, unique = T, method = "bray", binary = F, ...) {
  sample_id_1 <- sample_id_2 <- i <- j <- NULL
  # make "dist" object with beta values
  rel_abundance_matrix <- rel_abundance_matrix(ta, sample_name = sample_id)
  betas_dist <- vegdist(rel_abundance_matrix, method = method, binary = binary, ...)

  # save number of betas in betas_dist in shortcut variable
  n <- attr(betas_dist, "Size")

  # make tibble with beta values if we want only unique sample pairs
  if (unique) {
    betas <- expand.grid(i = 1:n, j = 1:n) %>%
      filter(i < j) %>%
      mutate(sample_id_1 = labels(betas_dist)[i]) %>%
      mutate(sample_id_2 = labels(betas_dist)[j]) %>%
      mutate(beta = betas_dist[n * (i - 1) - i * (i - 1) / 2 + j - i]) %>%
      select(-i, -j)

    # make tibble with beta values if we want all sample pairs (redundant!)
  } else {
    betas <- as.matrix(betas_dist) %>%
      as_tibble() %>%
      mutate(sample_id_1 = attr(betas_dist, "Labels")) %>%
      gather(key = "sample_id_2", value = "beta", -sample_id_1)
  }

  # add sample info to betas table
  betas <- ta$samples %>%
    `names<-`(names(.) %>% str_c("_2")) %>%
    right_join(betas, by = "sample_id_2", multiple = "all")
  betas <- ta$samples %>%
    `names<-`(names(.) %>% str_c("_1")) %>%
    right_join(betas, by = "sample_id_1", multiple = "all")

  # return betas table
  betas
}

#' Get prevalences of taxa in general or per condition
#'
#' Returns a tidy table of prevalences: taxon presence counts in samples,
#' overall or per condition.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#'
#' @param ta A tidytacos object.
#' @param condition A string denoting a categorical variable in the sample table.
#' @param pres_abs Whether to resort to presence/absence screening.
#' @export
prevalences <- function(ta, condition = NULL, pres_abs = F) {
  abundances_extended <-
    ta$counts %>%
    filter(count > 0) %>%
    left_join(ta$samples, by = "sample_id")

  if (is.null(condition)) {
    abundances_extended %>%
      count(taxon_id) %>%
      rename(prevalence = n)
  } else if (pres_abs) {
    condition <- sym(condition)

    abundances_extended %>%
      select(taxon_id, sample_id, !!condition) %>%
      mutate(presence = "present") %>%
      complete(nesting(!!condition, sample_id), taxon_id, fill = list(presence = "absent")) %>%
      count(taxon_id, !!condition, presence) %>%
      complete(taxon_id, !!condition, presence, fill = list(n = 0))
  } else {
    condition <- sym(condition)

    abundances_extended %>%
      count(taxon_id, !!condition) %>%
      rename(prevalence = n) %>%
      complete(taxon_id, !!condition, fill = list(prevalence = 0))
  }
}

#' Get mean relative abundances of taxa in general or per condition
#'
#' Returns tidy table with average relatively abundances of taxa, overall or per
#' condition.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#'
#' @param ta A tidytacos object.
#' @param condition A string representing a categorical variable to compute the relative abundances in every option of the variable
#' @export
mean_rel_abundances <- function(ta, condition = NULL) {
  # if rel_abundance not present: add and remove on exit
  if (!"rel_abundance" %in% names(ta$counts)) {
    ta <- add_rel_abundance(ta)
  }

  if (is.null(condition)) {
    ta$counts %>%
      select(sample_id, taxon_id, rel_abundance) %>%
      complete(sample_id, taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id) %>%
      summarize(mean_rel_abundance = mean(rel_abundance)) %>%
      ungroup()
  } else {
    condition <- sym(condition)

    ta$counts %>%
      left_join(ta$samples, by = "sample_id") %>%
      select(!!condition, sample_id, taxon_id, rel_abundance) %>%
      complete(nesting(!!condition, sample_id), taxon_id, fill = list(rel_abundance = 0)) %>%
      group_by(taxon_id, !!condition) %>%
      summarize(mean_rel_abundance = mean(rel_abundance)) %>%
      ungroup()
  }
}

#' Get all data in one single table
#'
#' `everything()` merges all three tidytacos tables into one very large table.
#'
#' @param ta A tidytacos object.
#' @export
everything <- function(ta) {
  # make and return large table
  ta$counts %>%
    left_join(ta$samples, by = "sample_id") %>%
    left_join(ta$taxa, by = "taxon_id")
}

#' Extract the sample table
#'
#' @param ta A tidytacos object.
#' @export
samples <- function(ta) ta$samples

#' Extract the taxon table
#'
#' @param ta A tidytacos object.
#' @export
taxa <- function(ta) ta$taxa

#' Extract the count table
#'
#' @param ta A tidytacos object.
#' @export
counts <- function(ta) ta$counts

#' Perform an adonis test
#'
#' This function executes the [adonis2][vegan::adonis2] function of the vegan package
#' and returns the result.
#'
#' Samples where one or more predictors are NA are removed.
#'
#' @importFrom stats as.formula
#' @param ta A tidytacos object.
#' @param predictors A character vector with predictors to include in the model.
#' @param permutations The number of permutations (more permutations takes
#'   longer but gives a more accurate p-value).
#' @inheritDotParams vegan::adonis2
#' @return An object of class "adonis" (see [adonis][vegan::adonis]).
#' @examples
#' res <- urt %>%
#'   perform_adonis(c("plate", "method"), by = "terms")
#' res
#' @export
perform_adonis <- function(ta, predictors, permutations = 999, ...) {
  counts_matrix <- ta %>%
    purrr::modify_at("samples", drop_na, one_of(predictors)) %>%
    process_sample_selection() %>%
    add_rel_abundance() %>%
    counts() %>%
    counts_matrix(value = "rel_abundance")

  formula_RHS <- paste0(predictors, collapse = " + ")

  metadata <- tibble(sample_id = rownames(counts_matrix)) %>%
    left_join(ta$samples, by = "sample_id")

  adonis2(
    as.formula(paste("counts_matrix", formula_RHS, sep = " ~ ")),
    metadata,
    permutations = permutations, ...
  )
}

#' Return a counts matrix
#'
#' This function returns a matrix with taxon counts; the rows are samples and
#' the columns are taxa.
#'
#' @param ta A tidytacos object.
#' @param sample_name The name of the variable in the sample table to use as row
#'   names (unquoted).
#' @param taxon_name The name of the variable in the taxon table to use as
#'   column names (unquoted).
#' @param value The name of the variable in the counts table to use as count
#' @param keep_empty_samples Should empty samples be included in the matrix?
#' @return A matrix with count values.
#'
#' @export
counts_matrix <- function(ta, sample_name = sample_id, taxon_name = taxon_id, value = count, keep_empty_samples = FALSE) {
  value <- rlang::enquo(value)
  if ("tidytacos" %in% class(ta)) {
    sample_name <- rlang::enquo(sample_name)
    taxon_name <- rlang::enquo(taxon_name)

    # create sample name if it doesn't exist
    if (!rlang::quo_name(sample_name) %in% names(ta$samples)) {
      ta$samples$sample_name <- ta$samples$sample_id
    }

    # create taxon name if it doesn't exist
    if (!rlang::quo_name(taxon_name) %in% names(ta$taxa)) {
      ta$taxa$taxon_name <- ta$taxa$taxon_id
    }

    tidy_count <- ta %>%
      change_id_samples(sample_id_new = {{ sample_name }}) %>%
      change_id_taxa(taxon_id_new = {{ taxon_name }}) %>%
      counts()
  } else {
    tidy_count <- ta
  }
  M <- tidy_count %>% tidy_count_to_matrix(value = {{ value }})
  # add empty samples
  if (
    keep_empty_samples &&
    any(ta %>% 
    add_total_count() %>% 
    samples() %>% 
    pull(total_count) == 0)) {
    empty_sample_ids <- ta %>%
      add_total_count() %>%
      samples() %>%
      filter(total_count == 0) %>%
      pull(sample_id)
    if (ncol(M) == 0) {
      stop("No taxa found in the counts table.")
    }
    empty_samples <- matrix(0, nrow = length(empty_sample_ids), ncol = ncol(M))
    rownames(empty_samples) <- empty_sample_ids
    colnames(empty_samples) <- colnames(M)
    M <- rbind(M, empty_samples)
  }
  M
}

#' Return a relative abundance matrix
#'
#' This function returns a relative abundance matrix; the rows are samples and
#' the column are taxa.
#'
#' @param ta A tidytacos object.
#' @param sample_name The name of the variable in the sample table to use as row
#'   names (unquoted).
#' @param taxon_name The name of the variable in the taxon table to use as
#'   column names (unquoted).
#' @return A matrix with abundance values.
#'
#' @export
rel_abundance_matrix <- function(ta, sample_name = sample_id, taxon_name = taxon_id) {
  if (
    !"tidytacos" %in% class(ta)
  ) {
    stop("first argument should be a tidytacos object")
  }

  sample_name <- rlang::enquo(sample_name)
  taxon_name <- rlang::enquo(taxon_name)

  # add relative abundances if not present
  if (!"rel_abundance" %in% names(ta$counts)) {
    ta <- add_rel_abundance(ta)
  }

  ta %>%
    counts_matrix(
      sample_name = {{ sample_name }}, taxon_name = {{ taxon_name }}, value = rel_abundance
    )
}

#' Return a list of taxon IDs per condition
#'
#' This function returns a named list of taxon_ids per distinct value of
#' a categorical column of the samples table.
#'
#' @param ta A tidytacos object.
#' @param condition The name of a variable in the sample table that contains a
#'   categorical value.
#' @param read_treshold The minimum read count to consider a taxon.
#'
#' @return A list of taxon_id vectors.
#'
#' @export
taxonlist_per_condition <- function(ta, condition, read_treshold = 0) {
  condition <- rlang::enquo(condition)
  condition_str <- rlang::quo_name(condition)

  # Allows input to be symbol or string
  if (!rlang::quo_is_symbol(condition)) {
    condition <- sym(condition_str)
  }

  if (!condition_str %in% names(ta$samples)) {
    error_message <-
      paste("Condition", condition_str, "not found in sample table.")
    stop(error_message)
  }

  distinct_conditions <- unique(ta$samples %>% pull(!!condition))

  select_taxa_for_condition <- function(var) {
    ta %>%
      filter_samples(!!condition == var) %>%
      filter_counts(count >= read_treshold)
  }
  ta_per_condition <- lapply(distinct_conditions, select_taxa_for_condition)
  names(ta_per_condition) <- distinct_conditions

  tt_all <- lapply(ta_per_condition, counts)

  lapply(tt_all, function(x) unique(x$taxon_id))
}
