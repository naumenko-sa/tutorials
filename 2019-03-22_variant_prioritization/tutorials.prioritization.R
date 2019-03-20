###############################################################################
##   Tutorial for CCM bioinformatics series, The Hospital for Sick Children, 2019-03-20
###############################################################################

# installation
# install.packages("tidyverse")

setwd("~/Desktop/teaching_n_learning/tutorials/2019-03-22_variant_prioritization/")

library(DBI)
library(tidyverse)

# download example database and additional files from:
# https://drive.google.com/drive/folders/0B_bLL10GwDnsVWl6cUhQZVRUTVk

dbname <- "gemini.db"
con <- dbConnect(RSQLite::SQLite(), dbname = dbname)
dbListTables(con)
dbListFields(con, "variants")

qrySample <- "select name from samples"

samples <- dbGetQuery(con, qrySample)
sample <- samples[[1]]

qryReport <- "select 
        v.ref as Ref,
        v.alt as Alt,
        v.chrom as Chrom,
        v.start+1  as Pos,
        v.impact as Variation,
        v.depth as Depth,
        v.qual_depth as Qual_depth,
        v.gene as Gene,
        g.ensembl_gene_id as Ensembl_gene_id,
        v.clinvar_disease_name as Clinvar,
        v.transcript as Ensembl_transcript_id,
        v.aa_length as AA_position,
        v.exon as Exon,
        v.pfam_domain as Pfam_domain,
        v.rs_ids as rsIDs,
        v.aaf_1kg_all as Maf_1000g,
        v.aaf_exac_all as Exac_maf,
        v.max_aaf_all as Maf_all,
        v.exac_num_het as Exac_het,
        v.exac_num_hom_alt as Exac_hom_alt,
        v.sift_score as Sift_score,
        v.polyphen_score as Polyphen_score,
        v.cadd_scaled as Cadd_score,
        v.aa_change as AA_change,
        from variants v, gene_detailed g
        where v.transcript=g.transcript and v.gene=g.gene";

variants <- as_tibble(dbGetQuery(con, qryReport))
            
variants_rare_missense <- filter(variants, Variation == "missense_variant",Maf_all < 0.01)
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
