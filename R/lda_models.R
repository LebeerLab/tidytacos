#' LDA model estimation
#'
#' `perform_lda()` estimates a Latent Dirichlet Allocation (LDA) model
#' on the counts matrix of a tidytacos object. The function returns
#' the estimated topics and terms, as well as the perplexity of the model.
#'
#' @param ta A tidytacos object.
#' @param k The number of topics to estimate.
#' @param min_prevalence The lowest percentage (0-1) of samples
#' taxa need to be present in for the taxa to be used in model estimation.
#' @param taxon The column name in the taxa table with taxa identifiers.
#' @param sample The column name in the sample table with sample identifiers.
#' @inheritDotParams topicmodels::LDA
#' @return A list of estimated topics, terms and the perplexity of the model.
#' @export
perform_lda <- function(ta, k,
                        min_prevalence = .05,
                        taxon = taxon_id, sample = sample_id, ...) {

  force_optional_dependency("topicmodels")
  taxon <- enquo(taxon)
  sample <- enquo(sample)

  matrix <- ta %>%
    add_prevalence(relative = TRUE) %>%
    filter_taxa(prevalence >= min_prevalence) %>%
    counts_matrix(sample = !!sample, taxon = !!taxon)

  model <- topicmodels::LDA(matrix, k, ...)
  lda <- topicmodels::posterior(model, matrix)

  results <- list(
    terms = lda$terms,
    topics = lda$topics,
    perplexity = topicmodels::perplexity(model)
  )

  results
}

#' Calculate LDA model perplexities for a range of topic numbers
#' 
#' `calculate_lda_perplexities()` estimates LDA models for a range of topic numbers 
#' and returns the perplexity of each model.
#' 
#' @param ta A tidytacos object.
#' @param k_range A range of topic numbers to estimate.
#' @param seed A seed for reproducibility.
#' @inheritDotParams perform_lda
#' @return A tibble with the perplexity of each model and the number of topics.
#' @export
calculate_lda_perplexities <- function(ta, k_range = 2:10, seed = 42, ...) {

  set.seed(seed)
  perplexity <- NULL

  perplexities <- purrr::map_dbl(k_range, function(.x, ...){
    perform_lda(ta, .x, ...) %>% `$`(perplexity)
  }, ...)

  tibble::tibble(k = k_range, perplexity = perplexities)
}

#' Aligns LDA topics across models
#' 
#' `align_lda_topics()` estimates LDA models for a range of topic numbers
#' and aligns the topics across models using the product and transport method.
#' 
#' @param ta A tidytacos object.
#' @param k_range A range of topic numbers to estimate.
#' @param seed A seed for reproducibility.
#' @param min_prevalence The lowest percentage (0-1) of samples
#' taxa need to be present in for the taxa to be used in model estimation.
#' @param taxon The column name in the taxa table with taxa identifiers.
#' @param sample The column name in the sample table with sample identifiers.
#' @inheritDotParams alto::run_lda_models
#' @return A list of estimated models, aligned topics using the product method
#' and aligned topics using the transport method.
#' @export
align_lda_topics <- function(ta, k_range = 2:10, seed = 42,
                             min_prevalence = .05,
                             taxon = taxon_id, sample = sample_id,
                             ...) {

  force_optional_dependency("alto", "devtools::install_github('lasy/alto')")
  set.seed(seed)

  lda_pars <- list()
  taxon <- dplyr::enquo(taxon)
  sample <- dplyr::enquo(sample)
  for (k in k_range) {
    lda_pars[[paste0("k", k)]] <- list(k = k)
  }

  matrix <- ta %>%
    add_prevalence(relative = TRUE) %>%
    filter_taxa(prevalence >= min_prevalence) %>%
    counts_matrix(sample = !!sample, taxon = !!taxon)

  lda_models <- alto::run_lda_models(
    data = matrix, lda_varying_params_lists = lda_pars,
    ...
  )

  product <- alto::align_topics(
    lda_models, method = "product"
  )

  transport <- alto::align_topics(
    lda_models, method = "transport"
  )
  list(models = lda_models, product = product, transport = transport)
}

#' Plot LDA topic alignment results
#'
#' `ldaplot_alignment()` plots the alignment of LDA topics across models.
#'
#' @param lda_alignment The output of `align_lda_topics()`.
#' @param show_product Whether to show the product or transport alignment.
#' @inheritDotParams alto::plot_alignment
#' @return A ggplot object.
#' @export
ldaplot_alignment <- function(lda_alignment, show_product = TRUE, ...) {

  res <- lda_alignment$transport

  if (show_product) {
    res <- lda_alignment$product
  }

  alto::plot(res, ...)
}

#' Plot LDA beta values
#'
#' `ldaplot_beta()` plots the beta values of LDA topics.
#'
#' @param lda_alignment The output of `align_lda_topics()`.
#' @param show_product Whether to show the product or transport alignment.
#' @inheritDotParams alto::plot_beta
#' @return A ggplot object.
#' @export
ldaplot_beta <- function(lda_alignment, show_product = TRUE, ...) {

  res <- lda_alignment$transport

  if (show_product) {
    res <- lda_alignment$product
  }

  alto::plot_beta(res, ...)

}