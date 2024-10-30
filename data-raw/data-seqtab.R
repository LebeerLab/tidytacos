library(dada2)

dada_ext <- system.file("extdata", package="dada2")

Fw <- paste0(dada_ext, c("/sam1F.fastq.gz", "/sam2F.fastq.gz"))
Rv <- paste0(dada_ext, c("/sam1R.fastq.gz", "/sam2R.fastq.gz"))

tax <- system.file("extdata", "example_train_set.fa.gz", package="dada2")
species_tax <- system.file("extdata", "example_species_assignment.fa.gz", package="dada2")


filtFw <- c(tempfile("filtFw1"), tempfile("filtFw2"))
filtRv <- c(tempfile("filtRv"), tempfile("filtRv2"))


out <- filterAndTrim(Fw,filtFw, Rv, filtRv)

errF <- learnErrors(filtFw, multithread=TRUE)
errR <- learnErrors(filtRv, multithread=TRUE)

dadaFw <- dada(filtFw, err=errF, multithread=TRUE)
dadaRv <- dada(filtRv, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFw, filtFw, dadaRv, filtRv, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "inst/extdata/dada2/seqtab")
taxa <- assignTaxonomy(seqtab, tax)
taxa <- addSpecies(taxa, species_tax)
saveRDS(taxa, "inst/extdata/dada2/taxa")
