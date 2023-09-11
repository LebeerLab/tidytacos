#' Add logratios
#'
#' This function computes pairwise logratio values between all taxa and adds
#' these to the tidytacos object in the form of a table called logratios.
#'
#' If max_taxa is greater than the number of taxa, the taxa with the highest
#' occurrence will be selected.
#'
#' IMPORTANT: this function add pseudocounts of one to all abundances before
#' calculating the logratios.
#'
#' @param ta A tidytacos object
#' @param max_taxa The maximum number of taxa to use
#'
#' @return A tidytacos object with an extra table logratios
#'
#' @export
add_logratios <- function(ta, max_taxa = 30) {

  if (nrow(ta$taxa) > max_taxa) {

    ta <- ta %>% add_occurrences()

    ta$taxa <-
      ta$taxa %>%
      arrange(desc(occurrence)) %>%
      mutate(keep = F) %>%
      {.$keep[1:max_taxa] <- T; .}

    ta <- ta %>% filter_taxa(keep) %>% 
    # cleaning up the keep column as it is not used further
          select_taxa(-keep)

  }

  abundances_complete <-
    ta$counts %>%
    complete(sample_id, taxon_id, fill = list(count = 0))

  ta$logratios <-
    left_join(
      abundances_complete %>%
        select(sample_id, taxon_id, count),
      abundances_complete %>%
        select(sample_id, ref_taxon_id = taxon_id, ref_abundance = count),
      by = "sample_id",
      multiple = "all"
    ) %>%
    mutate(
      taxon_ids = str_c(taxon_id, ref_taxon_id, sep = "_"),
      logratio = log10((count + 1) / (ref_abundance + 1))
    ) %>%
    select(- count, - ref_abundance)

  ta

}

#' Perform compositional differential abundance analysis
#'
#' This function performs a differential abundance test for all pairwise ratios
#' between taxa.
#'
#' A table called taxon_pairs will be added to the tidytacos object, with
#' for each pair of a taxon and a reference taxon, the differential abundance of
#' the taxon between the two conditions (with respect to the reference taxon).
#' The test that is performed is a Wilcoxon rank sum test and the test statistic
#' that is reported is the two-sample Hodges–Lehmann estimator (the median of
#' all pairwise differences between the samples).
#'
#' It is possible to supply the conditions to compare through the conditions
#' argument. Other conditions than the two supplied will be removed from the
#' data.
#'
#' @param ta A tidytacos object
#' @param condition A binary variable in the sample table (unquoted)
#' @param conditions A character vector with exactly two categories of the
#'   condition variable
#' @param max_taxa The maximum number of taxa to use
#'
#' @return A tidytacos object with an extra table taxon_pairs
#'
#' @export
add_codifab <- function(ta, condition, conditions = NULL, max_taxa = 30) {

  ta_sub <- ta

  condition <- rlang::enquo(condition)
  if (! rlang::f_text(condition) %in% names(ta$samples)) {
    stop("condition field does not exist in sample table")
  }

  ta_sub$samples <- ta_sub$samples %>% mutate(condition = !! condition)
  if (is.null(conditions)) {
    conditions <- unique(ta_sub$samples$condition)
  } else {
    conditions <- unique(conditions)
  }
  a_vs_b <- paste0(conditions[1], "_vs_", conditions[2])

  if (! length(conditions) == 2) {
    stop("there need to be exactly two conditions")
  }
  if (! all(conditions %in% unique(ta_sub$samples$condition))) {
    stop("one or both conditions not found")
  }

  ta_sub <- filter_samples(ta_sub, condition %in% conditions)

  # if logratios not present: add
  if (! "logratios" %in% names(ta_sub)) {
    ta_sub <- add_logratios(ta_sub, max_taxa = max_taxa)
  }

  ta$taxon_pairs <-
    ta_sub$logratios %>%
    filter(taxon_id != ref_taxon_id) %>%
    left_join(ta_sub$samples, by = "sample_id") %>%
    group_by(taxon_ids, taxon_id, ref_taxon_id) %>%
    summarize(
      wilcox = list(wilcox.test(
        x = logratio[condition == conditions[1]],
        y = logratio[condition == conditions[2]],
        conf.int = T, exact = F
      )),
      a_vs_b = purrr::map_dbl(wilcox, ~ .[["estimate"]]),
      wilcox_p = purrr::map_dbl(wilcox, ~ .[["p.value"]])
    ) %>%
    ungroup() %>%
    mutate(a_vs_b = 10 ^ a_vs_b) %>%
    rename(!! a_vs_b := a_vs_b)

  ta

}

#' Generate a compositional differential abundance plot
#'
#' This function returns a plot to visualize differential abundance of taxa
#' between conditions, compared to all other taxa as references. These
#' differential abundances should already have been calculated with
#' [add_codifab].
#'
#' Significance of tests is determined by capping the false discovery rate at
#' 10%, using the method of Benjamini and Yekutieli, which is developed for
#' non-independent tests. See [p.adjust].
#' @importFrom stats p.adjust median
#' @param ta A tidytacos object
#' @param diffabun_var The variable with differential abundances in the
#'   taxon_pair table
#'
#' @return A ggplot object
#'
#' @export
codifab_plot <- function(ta, diffabun_var) {

  if (! "taxon_name" %in% names(ta$taxa)) {
    ta <- add_taxon_name(ta)
  }

  diffabun_var <- rlang::enquo(diffabun_var)
  if (is.null(ta$taxon_pairs)) {
    stop("Please first run add_codifab() to generate the taxon pair comparisons.")
  } else if (! rlang::f_text(diffabun_var) %in% names(ta$taxon_pairs)) {
    stop(paste0(rlang::f_text(diffabun_var)," is not an existing comparison in the taxon_pairs table."))
  }

  taxon_pairs <-
    ta$taxon_pairs %>%
    mutate(wilcox_p = p.adjust(wilcox_p, "BY")) %>%
    left_join(
      ta$taxa %>% select(taxon_id, taxon = taxon_name), by = "taxon_id"
    ) %>%
    left_join(
      ta$taxa %>% select(ref_taxon_id = taxon_id, ref_taxon = taxon_name),
      by = "ref_taxon_id"
    ) %>%
    mutate(
      direction = if_else(!! diffabun_var > 1, "+", "-"),
      sign = wilcox_p < 0.10
    )

  taxa_ordered <-
    taxon_pairs %>%
    group_by(taxon) %>%
    summarize(median_diffabun = median(!! diffabun_var)) %>%
    arrange(median_diffabun) %>%
    pull(taxon)

  taxon_pairs %>%
    mutate_at(c("taxon", "ref_taxon"), factor, levels = taxa_ordered) %>%
    ggplot(aes(x = ref_taxon, y = taxon, fill = !! diffabun_var)) +
    geom_tile() +
    geom_text(
      aes(label = if_else(sign, direction, ""), col = direction), size = 2
    ) +
    scale_color_manual(values = c("+" = "black", "-" = "white"), guide = 'none') +
    scale_fill_continuous(trans = "log10") +
    xlab("reference taxon") +
    theme_minimal() +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
    )

}

#' Add compositional principal components to the sample table
#' @importFrom stats prcomp
#' @param ta a tidytacos object.
#' @export
add_copca <- function(ta) {

  # if logratios not present: add temporarily
  logratios_tmp <- ! "logratios" %in% names(ta)
  if (logratios_tmp) ta <- add_logratios(ta)

  logratio_matrix <-
    ta$logratios %>%
    select(taxon_ids, sample_id, logratio) %>%
    spread(key = taxon_ids, value = logratio) %>%
    {
      m <- as.matrix(.[, -1])
      row.names(m) <- .$sample_id
      m
    }

  pca <- prcomp(logratio_matrix[, colSums(logratio_matrix) != 0], scale. = T)
  samples_pca <- tibble(
    sample_id = rownames(pca$x),
    pca_1 = unname(pca$x[, 1]),
    pca_2 = unname(pca$x[, 2])
  )

  # add PCA dimensions to sample table
  ta$samples <- ta$samples %>% left_join(samples_pca, by = "sample_id")

  # cleanup
  if (logratios_tmp) ta$logratios <- NULL

  ta

}
