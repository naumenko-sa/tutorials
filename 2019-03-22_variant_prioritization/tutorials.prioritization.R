###############################################################################
##   Tutorial for CCM bioinformatics series, The Hospital for Sick Children, 2019-03-20
###############################################################################

# installation
# install.packages("tidyverse")
# install.packages("googledrive")

library(tidyverse)

setwd("~/Desktop/teaching_n_learning/2019-03-22_variant_prioritization/")

# download ashkenazim.small_variants.csv
# https://drive.google.com/open?id=1AiiITUwjlYQIssTWRZ-DTVS3a0_PFKJO 

library(googledrive)
drive_find(n_max=50, type = "csv")
#drive_download("~/tutorials/2019_03_variant_prioritization/ashkenazim.small_variants.csv")
id <- "1AiiITUwjlYQIssTWRZ-DTVS3a0_PFKJO"
drive_download(as_id(id), overwrite = T)
drive_download(as_id("https://drive.google.com/open?id=1AiiITUwjlYQIssTWRZ-DTVS3a0_PFKJO"), overwrite = T)

variants <- read_csv("ashkenazim.small_variants.csv")
            
variants_rare_missense <- filter(variants, Variation == "missense_variant", Maf_all < 0.01)

variants_potentially_deleterious <- filter(variants_rare_missense, Polyphen_score > 0.9)

#explore variation
effects <- count(variants, "Variation")
op <- par(mar = c(15,4,4,2) + 0.1)
barplot(effects$freq, names.arg = effects$Variation, las=2)
par(op)

#remove intronic variants
effects_no_intronic <- effects[effects$Variation != 'intron_variant',]
op <- par(mar = c(15,4,4,2) + 0.1)
barplot(effects_no_intronic$freq, names.arg = effects_no_intronic$Variation, las = 2)
par(op)

qryGene_detailed <- "select * from gene_detailed"
gene_detailed <- dbGetQuery(con, qryGene_detailed)

qryGene_summary <- "select * from gene_summary"
gene_summary <- dbGetQuery(con, qryGene_summary)

# additional information
# gene descriptions
gene_descriptions <- read_tsv("ensembl_w_description.txt")
variants <- left_join(variants, gene_descriptions, by = c("Ensembl_gene_id" = "ensembl_gene_id"))

# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/README_fordist_cleaned_exac_r03_z_data_pLI_2016_01_13.txt
exac_scores <- read_tsv("exac_scores.txt")
variants <- left_join(variants, exac_scores, by = c("external_gene_name" = "gene"))

qryAllTranscriptsAffected <- "select v.chrom, 
                                v.start+1,
                                v.end, 
                                v.ref,
                                v.alt,
                                v.gene,
                                vi.gene,
                                v.transcript,
                                vi.transcript,
                                vi.aa_change,
                                vi.aa_length 
                            from 
                                variant_impacts vi, 
                                variants v 
                            where 
                                vi.variant_id = v.variant_id and
                                v.is_coding = 1"
all_coding_effects = dbGetQuery(con,qryAllTranscriptsAffected)

dbDisconnect(con)
