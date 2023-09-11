# Prepare tidytacos object for visualization by barplot.
#
# Clusters samples, adds color groups and relative abundances.
#
# @param ta a tidytacos object
# @param n an integer
#
prepare_for_bp <- function(ta, n = 12, extended = TRUE) {

  # add sample_clustered if not present
  if (!"sample_clustered" %in% names(ta$samples)) {
    ta <- add_sample_clustered(ta)
  }

  # add taxon_name_color if not present
  if (!"taxon_name_color" %in% names(ta$taxa)) {
    ta <- add_taxon_name_color(ta, n = n)
  }

  # add relative abundances if not present
  if (!"rel_abundance" %in% names(ta$counts)) {
    ta <- add_rel_abundance(ta)
  }

  # optional extension (not used by sample bp)
  if (extended) {
    ta <- ta %>% everything()
  }
  ta
}

#' Return a bar plot of the samples
#'
#' @export
tacoplot_stack <- function(ta, n = 12, x = sample_clustered, geom_bar = T) {
  # convert promise to formula
  x <- enquo(x)

  warning_message_label = paste0("Label \'", quo_name(x),"\' not found in the samples table.")
  warning_message_aggregate = "Sample labels not unique, samples are aggregated."
  if (quo_name(x) != "sample_clustered" &&
    !is.element(quo_name(x), names(ta$samples))
  ) {
    # Warning, so tidy functions can be performed on the label
    warning(warning_message_label)
  }

  if (quo_name(x) != "sample_clustered" &&
    length(unique(ta$samples %>% pull(!!x))) < nrow(ta$samples)
  ) {
    warning(warning_message_aggregate)
  }

  # make plot and return
  plot <- prepare_for_bp(ta, n) %>%
    ggplot(aes(
      x = forcats::fct_reorder(!!x, as.integer(sample_clustered)),
      y = rel_abundance, fill = taxon_name_color)) +
    scale_fill_brewer(palette = "Paired", name = "Taxon") +
    xlab("sample") +
    ylab("relative abundance") +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white")
    )

  # add geom_bar if requested
  if (geom_bar) {
    plot <- plot + geom_bar(stat = "identity")
  }

  plot
}

#' Return an interactive bar plot of the samples
#'
#' @param ta A tidytacos object.
#' @param n An integer, representing the amount of colors used to depict
#'   different taxa.
#' @param x A string, representing the column name used to label and cluster the
#'   samples on.
#'
#' @export
tacoplot_stack_ly <- function(ta, n = 12, x = sample_clustered) {
  force_optional_dependency("plotly")
  # convert promise to formula
  x <- enquo(x)

  # wrap in eval and quosure shenannigans
  plot <- rlang::eval_tidy(rlang::quo_squash(
    quo({
      # make plot and return
      prepare_for_bp(ta, n) %>%
        plotly::plot_ly(
          x = ~forcats::fct_reorder(!!x, as.integer(sample_clustered)),
          y = ~rel_abundance,
          color = ~taxon_name_color,
          colors = palette_xgfs,
          name = ~taxon_name_color,
          hovertemplate = paste(
            "<b>%{x}</b>",
            "<br>%{y:.2%}<br>"
          ),
          type = "bar"
        ) %>%
        plotly::layout(
          barmode = "stack",
          xaxis = list(title="Sample"),
          yaxis = list(title="Relative Abundance")
        )
    })
  ))
  plot
}

#' Return an interactive pcoa plot of the samples
#'
#' @param ta A tidytacos object.
#' @param x A string, representing the column name used to color the sample
#'   groups on.
#' @param palette A vector of colors, used as the palette for coloring sample
#'   groups.
#'
#' @export
tacoplot_ord_ly <- function(ta, x=NULL, samplenames = sample_id, palette = NULL, title = "PCOA plot") {
  force_optional_dependency("plotly")
  # convert promise to formula
  if (is.null(x)) {
    stop("Argument x missing. Please supply the name of a categorical value, to be used as the color for the pcoa plot.")
  }
  
  x <- enquo(x)
  samplenames <- enquo(samplenames)

  # fallback to default palette
  if (is.null(palette)) {
    palette <- palette_paired
  }

  # prepare pcoa if needed
  if (!all(c("pcoa1", "pcoa2") %in% names(ta$samples))) {
    ta <- add_pcoa(ta)
  }

  plot <- rlang::eval_tidy(rlang::quo_squash(
    quo({
      ta$samples %>%
        plotly::plot_ly(
          x = ~pcoa1,
          y = ~pcoa2,
          color = ~!!x,
          colors = palette_paired,
          text = ~!!samplenames,
          hovertemplate = paste("<i>%{text}</i>"),
          type = "scatter",
          mode = "markers"
        ) %>%
        plotly::layout(
          title = title,
          yaxis = list(zeroline = F),
          xaxis = list(zeroline = F)
        )
    })
  ))

  plot
}

#' Return a pcoa plot of the samples
#'
#' @param ta A tidytacos object.
#' @param x A string, representing the column name used to color the sample
#'   groups on.
#' @param palette A vector of colors, used as the palette for coloring sample
#'   groups.
#'
#' @export
tacoplot_ord <- function(ta, x=sample_id, palette = NULL, title = "PCOA plot") {

  x <- enquo(x)
  
  error_message = paste0("Label \'", quo_name(x),"\' not found in the samples table.")
  if(!is.element(quo_name(x), names(ta$samples))) {
    stop(error_message)
  }

  if (quo_name(x) == "sample_id") {
    x <- NULL
  }
  
  # fallback to default palette
  if (is.null(palette)) {
    palette <- palette_paired
  }

  # prepare pcoa if needed
  if (!all(c("pcoa1", "pcoa2") %in% names(ta$samples))) {
    ta <- add_pcoa(ta)
  }

  ta$samples %>% ggplot(aes(x=pcoa1, y=pcoa2, color=!!x)) + 
    geom_point() + 
    theme_classic() +
    ggtitle(title)

}

#' Return a visualization designed for a small number of samples
#'
#' @export
tacoplot_zoom <- function(ta, sample = sample_id, n = 15, nrow = NULL) {
  ta <- prepare_for_bp(ta, n, extended = FALSE)

  sample <- rlang::enexpr(sample)
  if (sample != rlang::expr(sample_id)) {
    ta <- change_id_samples(ta, sample_id_new = !!sample)
  }

  data <-
    ta %>%
    everything() %>%
    group_by(sample_id) %>%
    arrange(desc(rel_abundance)) %>%
    slice(1:n) %>%
    ungroup() %>%
    arrange(sample_id, rel_abundance) %>%
    mutate(row = 1:n())

  data %>%
    ggplot(aes(x = row, y = rel_abundance, fill = taxon_name_color)) +
    geom_col() +
    facet_wrap(~sample_id, scales = "free", nrow) +
    coord_flip() +
    theme_bw() +
    scale_x_continuous(
      breaks = data$row,
      labels = data$taxon_name,
      expand = c(0, 0)
    ) +
    scale_fill_brewer(palette = "Paired", name = "taxon") +
    xlab("taxon name") +
    ylab("relative abundance")
}

#' Return a venn diagram of overlapping taxon_ids between conditions
#'
#' @param ta A tidytacos object.
#' @param condition The name of a variable in the samples table that contains a
#'   categorical value.
#'
#' @export
tacoplot_venn <- function(ta, condition, ...) {

  force_optional_dependency("ggVennDiagram")

  condition <- enquo(condition)
  ltpc <- taxonlist_per_condition(ta, !!condition)
  ggVennDiagram::ggVennDiagram(ltpc, ...)

}

palette_paired <- c(
  "#e8e8e8", # light grey
  "#a6cee3", # light blue
  "#1f78b4", # dark blue
  "#b2df8a", # light green
  "#33a02c", # dark green
  "#fb9a99", # light red
  "#e31a1c", # dark red
  "#fdbf6f", # light orange
  "#ff7f00", # dark orange
  "#cab2d6", # light purple
  "#6a3d9a", # dark purple
  "#ffff99", # light brown
  "#b15928" # dark brown
)

palette_xgfs <- c(
  # source: http://tsitsul.in/blog/coloropt/
  "#bdbdbd",
  "#00a76c",
  "#878500",
  "#00c6f8",
  "#5954d6",
  "#ff9287",
  "#b24502",
  "#d163e6",
  "#00bbad",
  "#006e00",
  "#008cf9",
  "#b80058",
  "#ebac23"
)
