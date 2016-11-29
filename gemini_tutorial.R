##########################################################
##   Tutorial for CCM bioinformatics series, Sick Kids
##########################################################

setwd("~/Desktop/teaching_n_learning/tutorials/")

library(RSQLite)
library(plyr)

dbname="NA12878-1-ensemble.db"

con = dbConnect(RSQLite::SQLite(),dbname=dbname)

dbListTables(con)

qrySample="select name from samples"

samples = dbGetQuery(con,qrySample)
sample = samples[[1]]

qryReport=paste0("select 
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
        v.aa_change as AA_change
        from variants v, gene_detailed g
        where v.transcript=g.transcript and v.gene=g.gene");

variants = dbGetQuery(con,qryReport)

variants.rare.missense = subset(variants,Variation == 'missense_variant' & Maf_all < 0.01)

variants.potentially_harmful = subset(variants.rare.missense, Polyphen_score > 0.9)

#explore variation
effects = count(variants,'Variation')
op <- par(mar = c(15,4,4,2) + 0.1)
barplot(effects$freq,names.arg = effects$Variation,las=2)
par(op)

#remove intronic variants
effects.no_intronic = effects[effects$Variation != 'intron_variant',]
op <- par(mar = c(15,4,4,2) + 0.1)
barplot(effects.no_intronic$freq,names.arg = effects.no_intronic$Variation,las=2)
par(op)

qryGene_detailed = "select * from gene_detailed"
gene_detailed = dbGetQuery(con,qryGenes)

qryGene_summary = "select * from gene_summary"
gene_summary = dbGetQuery(con,qryGene_summary)

#additional information
#gene descriptions
#reference_tables_path="~/Desktop/reference_tables"
reference_tables_path=getwd()

gene_descriptions = read.delim2(paste0(reference_tables_path,"/ensembl_w_description.txt"), stringsAsFactors=FALSE)
variants = merge(variants,gene_descriptions,by.x = "Ensembl_gene_id",by.y = "ensembl_gene_id",all.x=T)

#ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/README_fordist_cleaned_exac_r03_z_data_pLI_2016_01_13.txt
exac_scores = read.delim(paste0(reference_tables_path,"/exac_scores.txt"), stringsAsFactors=F)
variants = merge(variants,exac_scores,all.x=T)

dbDisconnect(con)
