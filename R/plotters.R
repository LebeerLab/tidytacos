# Prepare tidytacos object for visualization by barplot.
#
# Clusters samples, adds color groups and relative abundances.
#
# @param ta a tidytacos object
# @param n an integer
#
prepare_for_bp <- function(ta, n = 12, extended = TRUE, order_by = NULL, aggregate=TRUE) {
  # custom order
  if (!is.null(order_by)) {
    ta$samples <- ta$samples %>% arrange(order_by)
    ta$samples$sample_clustered <- as.factor(ta$samples[[order_by]])
  }


  # add taxon_name_color if not present
  if (!"taxon_name_color" %in% names(ta$taxa)) {
    ta <- add_taxon_name_color(ta, n = n)
  }

  # add relative abundances if not present
  if (!"rel_abundance" %in% names(ta$counts)) {
    ta <- add_rel_abundance(ta)
  }

if (aggregate){
      # aggregate taxa by taxon_name_color
      # to prevent many streaks for the 'Other' taxa
      ta <- ta %>%
      set_rank_names("taxon_name_color") %>%
      aggregate_taxa(rank="taxon_name_color")
    } 
 
  # add sample_clustered if not present
  if (!"sample_clustered" %in% names(ta$samples)) {
    ta <- add_sample_clustered(ta)
  }

  # optional extension (not used by sample bp)
  if (extended) {
       ta <- ta %>% everything()
  }
  ta
}

# Calculate beta diversity statistics to display on an ordination plot.
#
#
# @param ta a tidytacos object
# @param x the variable of interest
# @param stat.method the statistic to calculate, choice from anosim or permanova
# @param distance the beta dissimilarity metric to use
#
get_ord_stat <- function(ta, x, stat.method, distance = distance) {
  stat.method <- tolower(stat.method)
  if (stat.method == "anosim") {
    stat.result <- perform_anosim(ta, !!x, distance = distance)
  } else if (stat.method == "permanova" || stat.method == "adonis") {
    res <- perform_adonis(ta, c(rlang::as_name(x)), distance = distance)
    stat.result <- list(statistic = res$R2[[1]], signif = res$`Pr(>F)`[[1]])
  } else {
    stop("stat.method not recognized. Please choose from anosim or permanova.")
  }
  stat.result
}

#' Return a bar plot of the samples
#'
#' Plots a stacked bar plot of the samples in the tidytacos object to inspect the taxonomic profile.
#'
#' @param ta A tidytacos object.
#' @param x The name of the column name used to represent samples on the x-axis
#' @param n An integer, representing the amount of colors used to depict
#' @param pie A boolean, whether or not to represent the profile in a pie chart.
#' Default is FALSE, as pie chart representations can be misleading to interpret.
#' @param order_by an optional column name to order the samples by.
#' For examples order_by=sample would order the x-axis by the sample names instead of by similar profiles.
#' @export
tacoplot_stack <- function(ta, n = 12, x = sample_clustered, pie = FALSE, order_by = NULL) {
  # convert promise to formula


  try(
    {
      if (is.function(x)) {
        x <- as.character(substitute(x))
      }
      x <- rlang::sym(x)
    },
    silent = TRUE
  ) # in case the column is given as a string
  x <- rlang::enquo(x)

  error_message_label <- paste0("Label \'", rlang::quo_name(x), "\' not found in the samples table.")
  error_message_pie <- "This visualization type is meant to be used for a single sample."
  warning_message_aggregate <- "Sample labels not unique, samples are aggregated."

  if (rlang::quo_name(x) != "sample_clustered" &&
    !is.element(rlang::quo_name(x), names(ta$samples))
  ) {
    # Warning, so tidy functions can be performed on the label
    stop(error_message_label)
  }

  if (quo_name(x) != "sample_clustered" &&
    length(unique(ta$samples %>% pull(!!x))) < nrow(ta$samples)
  ) {
    warning(warning_message_aggregate)
  }

  if (pie & length(ta$samples) > 1) {
    stop(error_message_pie)
  }

  # make plot and return
  plot <- prepare_for_bp(ta, n, extended = T, order_by = order_by) %>%
    ggplot(aes(
      x = forcats::fct_reorder(!!x, as.integer(sample_clustered)),
      y = rel_abundance, fill = taxon_name_color
    )) +
    list(
      geom_bar(stat = "identity"),
      scale_fill_brewer(palette = "Paired", name = "Taxon"),
      if (pie) coord_polar("y", start = 0, clip = "off"),
      xlab("sample"),
      ylab("relative abundance"),
      theme(
        axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white")
      ),
      if (pie) xlab(""),
      if (pie) theme(axis.text.x = element_blank())
    )

  # Add > 12 colors if asked for
  if (n > 12) {
    force_optional_dependency("RColorBrewer")
    suppressMessages(
      plot <- plot +
        scale_fill_manual(values = colorRampPalette(palette_xgfs)(n))
    )
  }

  plot
}

#' Return an interactive bar plot of the samples
#'
#' Plots an interactive stacked bar plot of the samples in the tidytacos object to inspect the taxonomic profile.
#'
#' @param ta A tidytacos object.
#' @param n An integer, representing the amount of colors used to depict
#'   different taxa.
#' @param x A string, representing the column name used to label the x-axis
#' @param order_by an optional column name to order the samples by.
#' For examples order_by=sample would order the x-axis by the sample names instead of by similar profiles.
#'
#' @export
tacoplot_stack_ly <- function(ta, n = 12, x = sample_clustered, order_by = NULL) {
  force_optional_dependency("plotly")
  # convert promise to formula
  x <- rlang::enquo(x)

  # wrap in eval and quosure shenannigans
  plot <- rlang::eval_tidy(rlang::quo_squash(
    quo({
      # make plot and return
      prepare_for_bp(ta, n, extended = T, order_by = order_by, aggregate = FALSE) %>%
        plotly::plot_ly(
          x = ~ forcats::fct_reorder(!!x, as.integer(sample_clustered)),
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
          xaxis = list(title = "Sample"),
          yaxis = list(title = "Relative Abundance")
        )
    })
  ))
  plot
}

#' Return an interactive ordination plot of the samples
#'
#' Creates an interactive ordination plot of the beta diversity of the samples in the tidytacos object.
#' This can be used to gauge the similarity between samples.
#'
#' @param ta A tidytacos object.
#' @param x A string, representing the column name used to color the sample
#'   groups on.
#' @param samplenames the column in the sample table with the samplenames, defaults to sample_id.
#' @param ord the ordination technique to use. Choice from pcoa, tsne and umap.
#' @param distance the distance algorithm to use, see [vegan::vegdist()].
#' @param dims the amount of dimensions to plot, 2 or 3.
#' @param stat.method the statistic to print on the figure, choice from mantel and anosim.
#' @param palette A vector of colors, used as the palette for coloring sample
#' @param title a string to display as title of the plot.
#'   groups.
#' @param ... Extra arguments to pass to the add_ord function.
#'
#' @export
tacoplot_ord_ly <- function(ta, x = NULL, samplenames = sample_id, ord = "pcoa", dims = 2,
                            distance = "bray", stat.method = NULL, palette = NULL, title = NULL, ...) {
  force_optional_dependency("plotly")


  if (is.null(title)) {
    title <- paste(toupper(ord), "plot")
  }

  # convert promise to formula
  x <- rlang::enquo(x)
  if (rlang::quo_is_null(x)) {
    stop("Argument x missing. Please supply the name of a categorical value, to be used as the color for the pcoa plot.")
  }
  samplenames <- rlang::enquo(samplenames)
  ordnames <- c("ord1", "ord2")
  if (dims == 3) {
    ordnames <- c(ordnames, "ord3")
  }
  # fallback to default palette
  if (is.null(palette)) {
    palette <- palette_paired
  }

  # Check for empty samples
  ta <- ta %>% remove_empty_samples()

  # prepare ord if needed
  if (!all(ordnames %in% names(ta$samples))) {
    ta <- add_ord(ta, distance = distance, method = ord, dims = dims, ...)
  }

  if (dims == 2) {
    plot <- rlang::eval_tidy(rlang::quo_squash(
      quo({
        ta$samples %>%
          plotly::plot_ly(
            x = ~ord1,
            y = ~ord2,
            color = ~ !!x,
            colors = palette,
            text = ~ !!samplenames,
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
  } else {
    plot <- rlang::eval_tidy(rlang::quo_squash(
      quo({
        ta$samples %>%
          plotly::plot_ly(
            x = ~ord1,
            y = ~ord2,
            z = ~ord3,
            color = ~ !!x,
            colors = palette,
            text = ~ !!samplenames,
            hovertemplate = paste("<i>%{text}</i>")
          ) %>%
          plotly::add_markers() %>%
          plotly::layout(
            title = title,
            yaxis = list(zeroline = F),
            xaxis = list(zeroline = F)
          )
      })
    ))
  }

  if (is.null(stat.method)) {
    return(plot)
  }
  stat <- get_ord_stat(ta, x, stat.method, distance = distance)
  plot %>% plotly::add_annotations(
    x = 0.1,
    y = 1,
    xref = "paper",
    yref = "paper",
    text = paste0(toupper(stat.method), ":\nR= ", signif(stat$statistic, 3), "\nP= ", signif(stat$signif, 3)),
    showarrow = F
  )
}



#' Return an ordination plot of the samples
#'
#' Creates an ordination plot of the beta diversity of the samples in the tidytacos object.
#' This can be used to gauge the similarity between samples.
#'
#' @param ta A tidytacos object.
#' @param x The column name used to color the sample groups on.
#' @param ord the ordination technique to use. Choice from pcoa, tsne and umap.
#' @param distance the distance algorithm to use, see [vegan::vegdist()].
#' @param stat.method the statistic to print on the figure, choice from mantel and anosim.
#' @param palette A vector of colors, used as the palette for coloring sample
#'   groups.
#' @param title a string to display as title of the plot.
#' @param ... Extra arguments to pass to the add_ord function.
#'
#' @examples
#'
#' tacoplot_ord(urt, x = location)
#'
#' # set plate to character, to avoid it being treated as a continuous variable
#' urt <- urt %>% mutate_samples(plate = as.character(plate))
#' tacoplot_ord(urt, x = plate, ord = "umap", distance = "aitchison", stat.method = "permanova")
#'
#' @export
tacoplot_ord <- function(ta, x = NULL, palette = NULL, ord = "pcoa", distance = "bray", stat.method = NULL, title = NULL, ...) {
  # convert promise to formula
  x <- rlang::enquo(x)
  if (rlang::quo_is_null(x)) {
    stop("Argument x missing. Please supply the name of a categorical value, to be used as the color for the pcoa plot.")
  }

  if (is.null(title)) {
    title <- paste(ord, "plot")
  }
  error_message <- paste0("Label \'", quo_name(x), "\' not found in the samples table.")
  if (!is.element(quo_name(x), names(ta$samples))) {
    stop(error_message)
  }

  if ((quo_name(x) == "sample_id" ||
    quo_name(x) == "sample") &&
    (length(ta$samples$sample_id) > 20)) {
    warning(
      paste0(
        "Ignoring sample_id or sample as colouring variable for this many samples. ",
        "\nUse tacoplot_ord_ly(x=<group_var>, samplenames=sample_id) ",
        "if you wish to view the samples by name or id."
      )
    )
    x <- NULL
    stat.method <- NULL
  }

  # fallback to default palette
  if (is.null(palette)) {
    palette <- palette_paired
  }

  # Check for empty samples
  ta <- ta %>% remove_empty_samples()

  # prepare pcoa if needed
  if (!all(c("ord1", "ord2") %in% names(ta$samples))) {
    ta <- add_ord(ta, distance = distance, method = ord, ...)
  }

  plt <- ta$samples %>% ggplot(aes(x = ord1, y = ord2, color = !!x)) +
    geom_point() +
    theme_classic() +
    ggtitle(title) +
    labs(subtitle = paste("Distance measure: ", distance))

  # calculate stats
  if (is.null(stat.method)) {
    return(plt)
  }
  stat.result <- get_ord_stat(ta, x, stat.method, distance = distance)

  plt +
    annotate("text",
      x = min(ta$samples$ord1) + 0.05, y = max(ta$samples$ord2) - 0.05,
      label = paste0(
        toupper(stat.method),
        ":\nR= ", signif(stat.result$statistic, 3), "\nP= ", signif(stat.result$signif, 3)
      )
    )
}

#' Return a visualization designed for a small number of samples
#'
#' @param ta A tidytacos object.
#' @param sample A string, representing the unique sample name of interest
#' @param n An integer, representing the amount of colors used to depict taxa
#' @param nrow An integer, representing the amount of rows in the facet_wrap
#'
#' @export
tacoplot_zoom <- function(ta, sample = sample_id, n = 15, nrow = NULL) {
  ta <- prepare_for_bp(ta, n, extended = FALSE)


  sample <- rlang::enexpr(sample)
  if (sample != rlang::expr(sample_id)) {
    ta <- change_id_samples(ta, sample_id_new = !!sample)
  }

  if (ta$samples %>% nrow() > 1) {
    stop("This visualization is meant to be used for a single sample.")
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
      labels = data$taxon_name_color,
      expand = c(0, 0)
    ) +
    {
      if (n <= 12) scale_fill_brewer(palette = "Paired", name = "taxon")
    } +
    {
      if (n > 12) scale_fill_manual(values = colorRampPalette(palette_xgfs)(n))
    } +
    xlab("taxon name") +
    ylab("relative abundance") +
    theme(
      legend.position = "none"
    )
}

#' Return a venn diagram of overlapping taxon_ids between conditions
#'
#' @param ta A tidytacos object.
#' @param condition The name of a variable in the samples table that contains a
#'   categorical value.
#' @inheritDotParams ggVennDiagram::ggVennDiagram
#' @examples
#' tacoplot_venn(urt, location)
#' @export
tacoplot_venn <- function(ta, condition, ...) {
  force_optional_dependency("ggVennDiagram")


  condition <- rlang::enquo(condition)
  ltpc <- taxonlist_per_condition(ta, !!condition)

  if ("show_intersect" %in% names(list(...))) {
    if (!"taxon_name" %in% colnames(ta$taxa)) {
      ta <- ta %>% add_taxon_name()
    }
    match_taxon_name <- function(taxid) {
      ta$taxa[which(ta$taxa$taxon_id %in% taxid), ] %>%
        dplyr::pull(taxon_name)
    }
    ltpc <- lapply(ltpc, match_taxon_name)
  }

  ggVennDiagram::ggVennDiagram(ltpc, ...)
}

#' Return an interactive venn diagram of overlapping taxon_ids between conditions
#'
#' @param ta A tidytacos object.
#' @param condition The name of a variable in the samples table that contains a
#'   categorical value.
#' @inheritDotParams ggVennDiagram::ggVennDiagram
#'
#' @export
tacoplot_venn_ly <- function(ta, condition, ...) {
  condition <- rlang::enquo(condition)

  if (!"taxon_name" %in% names(ta$taxa)) {
    ta <- ta %>% add_taxon_name()
  }

  ta %>% tacoplot_venn(!!condition, show_intersect = TRUE, ...)
}

#' Return an euler diagram of overlapping taxon_ids between conditions
#'
#' @param ta A tidytacos object.
#' @param condition The name of a variable in the samples table that contains a
#'   categorical value.
#' @param shape shape to plot the groups in; choice from circle or ellipse
#' @inheritDotParams eulerr::euler
#' @export
tacoplot_euler <- function(ta, condition, shape = "ellipse", ...) {
  force_optional_dependency("eulerr")

  condition <- rlang::enquo(condition)
  ltpc <- taxonlist_per_condition(ta, !!condition)

  fit <- eulerr::euler(ltpc, shape)
  plot(fit, quantities = T, ...)
}

#' Return a boxplot of every alpha metric per group in the samples table of a tidytaco object.
#' If no alpha metrics are present, all available ones are added.
#'
#' @param ta A tidytacos object.
#' @param group_by The name of a variable in the samples table on which to group the samples.
#' @param compare_means Add the result of a statistical test to the plot, comparing the means of the groups. Default is FALSE.
#' @param keep_empty_samples Whether to discard samples not containing any counts or not. By default these are removed.
#' @param subsample_metric if the alpha diversities are precalculated with `add_subsampled_alpha()`, choose here "mean" or "median" to represent the alpha diversity.
#' @inheritDotParams ggpubr::stat_compare_means
#' @export
tacoplot_alphas <- function(ta, group_by, compare_means = FALSE, keep_empty_samples=FALSE, subsample_metric=NULL, ...) {
  value <- NULL

  group_by <- rlang::enquo(group_by)
  if (rlang::quo_is_missing(group_by)) {
    stop("Argument group_by missing. Please supply the name of a categorical value, to be used as the grouping variable or specify NULL if you just want to plot the alpha values of all samples.")
  }
  ta_tmp <- ta
  clean_alpha_metrics <- sapply(alpha_metrics, tolower, USE.NAMES = F)

  if (rlang::quo_is_null(group_by)) {
    group_by <- "all.samples"
    ta_tmp$samples$all.samples <- "all.samples"
  }

  if (!is.null(subsample_metric)){
    if (!subsample_metric %in% c("mean","median")){
      stop("'subsample_metric' can only be 'mean' or 'median'")
    }
    if (!any(colnames(ta_tmp$samples)%>% startsWith(paste0(subsample_metric,"_")))) {
      stop("Please first run add_alphas using subsample=TRUE")
    }
    ta_tmp$samples <- ta_tmp$samples %>% rename_with(~str_remove(., paste0(subsample_metric,"_"))) 
  }

  if (!any(clean_alpha_metrics %in% colnames(ta_tmp$samples))) {
    ta_tmp <- add_alphas(ta_tmp, keep_empty_samples=keep_empty_samples)
  }

  plt <- ta_tmp$samples %>%
    pivot_longer(any_of(clean_alpha_metrics)) %>%
    ggplot(aes(x = !!group_by, y = value, fill = !!group_by)) +
    geom_violin() +
    geom_jitter(alpha = 0.1) +
    facet_wrap(~name, scales = "free") +
    theme_classic() +
    theme(
      strip.background = element_rect(
        fill = alpha("lightblue", 0.4),
        linewidth = 0.5
      )
    )

  if (!compare_means) {
    return(plt)
  }

  force_optional_dependency("ggpubr")
  plt + ggpubr::stat_compare_means(...)
}

#' Return a heatmap of prevalence of taxa in groups of samples
#'
#' Return a heatmap of all taxa above a certain threshold prevalence
#' per condition, clusters them and compares prevalences with a fisher_test.
#' @examples
#' urt %>%
#'   aggregate_taxa(rank = "order") %>%
#'   tacoplot_prevalences(location, cutoff = .1,
#'   treeheight_row = 0, cutree_rows = 4,
#'   fontsize = 6, cellwidth = 15)
#'
#' @param ta A tidytacos object.
#' @param condition The row name of the condition which rrevalences are to be compared.
#' @param cutoff The minimum prevalence of a taxon to be included in the heatmap.
#' @param fisher Run a fisher test on the relative prevalences in each condition
#' and plot the resulting adjusted p-values as *(<.05), **(<.01), ***(<.001) or ****(<.0001).
#' @param adjp_method The method to adjust the p-values, see [rstatix::adjust_pvalue()].
#' @inheritDotParams pheatmap::pheatmap
#' @export
tacoplot_prevalences <- function(ta, condition, cutoff = 0.1, fisher = T, adjp_method = "fdr", ...) {
  force_optional_dependency("pheatmap")
  force_optional_dependency("rstatix")
  fisher_p <- NULL
  condition <- rlang::enquo(condition)

  prevalences <- ta %>%
    `if`(!"taxon_name" %in% ta$taxa, add_taxon_name(.), .) %>%
    add_prevalence(condition = rlang::as_name(condition), relative = T, fisher_test = fisher) %>%
    taxa()

  prevalences.M <- prevalences %>%
    dplyr::select(starts_with("prevalence_in_")) %>%
    as.matrix()
  if (fisher) {
    fp <- prevalences %>%
      dplyr::select(fisher_p) %>%
      rstatix::adjust_pvalue(method = adjp_method)
    fp_sig <- dplyr::case_when(
      fp < 0.0001 ~ " (****)",
      fp < 0.001 ~ " (***)",
      fp < 0.01 ~ " (**)",
      fp < 0.05 ~ " (*)",
      TRUE ~ ""
    )

    row.names(prevalences.M) <- str_c(prevalences$taxon_name, fp_sig)
  } else {
    row.names(prevalences.M) <- prevalences$taxon_name
  }

  prevalences.M <- prevalences.M[prevalences.M %>% rowSums() / ncol(prevalences.M) > cutoff, ]

  pheatmap::pheatmap(prevalences.M, ...)
}


#' Return a scree plot to visualize the eigenvalues of the PCA.
#'
#'
#' @param ta A tidytacos object.
#' @inheritDotParams factoextra::fviz_eig
#' @examples
#' urt %>% tacoplot_scree()
#' @export
tacoplot_scree <- function(ta, ...) {
  force_optional_dependency("factoextra")

  if (!"pca" %in% names(ta)) {
    ta <- suppressWarnings(add_copca(ta))
  }

  factoextra::fviz_eig(ta$pca, ...)
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
