test_that("Can extract topics from taco", {
    skip_if_not_installed("topicmodels")
    TOPICS <- 4
    
    urt.k4 <- urt %>% perform_lda(k=TOPICS, sample_id = sample)
    
    expect_equal(ncol(urt.k4$topics), TOPICS)
    expect_equal(nrow(urt.k4$terms), TOPICS)
})

test_that("Can align range of lda topics from taco", {

    skip_if_not_installed("alto")
    skip_if_not_installed("topicmodels")
    
    WARN_MESSAGE <- "Using default value 'VEM' for 'method' LDA parameter."
    TOPICS_MIN <- 2
    TOPICS_MAX <- 4

    urt_g <- urt %>% 
      aggregate_taxa(rank="genus") %>% 
      add_taxon_name()
    
    expect_message(
        expect_message(
            expect_message(
                lda_alignment <- align_lda_topics(
                    urt_g, k_range=TOPICS_MIN:TOPICS_MAX, taxon=taxon_name),
                WARN_MESSAGE
            ), WARN_MESSAGE
        ), WARN_MESSAGE
    )

    expect_equal(length(lda_alignment$models), TOPICS_MAX - TOPICS_MIN + 1)
    expect_equal(levels(lda_alignment$product@topics$m), c("k2", "k3", "k4"))
    expect_equal(levels(lda_alignment$transport@topics$m), c("k2", "k3", "k4"))

})

test_that("Can plot alignment and betas of lda topics", {

    skip_if_not_installed("alto")
    skip_if_not_installed("topicmodels")
    
    WARN_MESSAGE <- "Using default value 'VEM' for 'method' LDA parameter."
    TOPICS_MIN <- 1
    TOPICS_MAX <- 2

    urt_g <- urt %>% 
      aggregate_taxa(rank="genus") %>% 
      add_taxon_name()
    
    expect_message(
        expect_message(
            lda_alignment <- align_lda_topics(
                urt_g, k_range=TOPICS_MIN:TOPICS_MAX, taxon=taxon_name),
                WARN_MESSAGE
        ), WARN_MESSAGE
    )

    align_plot <- ldaplot_alignment(lda_alignment)
    align_plot_tr <- ldaplot_alignment(lda_alignment, show_product=FALSE)
    lda_betas <- ldaplot_beta(lda_alignment)

    vdiffr::expect_doppelganger("LDA alignment product", align_plot)
    vdiffr::expect_doppelganger("LDA alignment transport", align_plot_tr)
    vdiffr::expect_doppelganger("LDA betas", lda_betas)

})
