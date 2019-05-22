###############################################################################
##   Tutorial for CCM bioinformatics series, The Hospital for Sick Children, 2019-03-20
###############################################################################
# use select and then Ctrl + Enter to source a function
filter_variants <- function(variants, proband, mom, dad){
    variants_filtered <- filter(variants, Zygosity.Ashkenazim_HG002 == proband, 
                                 Zygosity.Ashkenazim_HG003 == mom,
                                 Zygosity.Ashkenazim_HG004 == dad)    
    return(variants_filtered)
}

# installation
# install.packages("tidyverse")
# install.packages("googledrive")

library(tidyverse)

setwd("~/Desktop/teaching_n_learning/2019-03-22_variant_prioritization/")

# download ashkenazim.small_variants.csv
# https://drive.google.com/open?id=15bYaktk50qD6qeL-LPgm5NSlk4-XFfYq

library(googledrive)
drive_find(n_max = 50, type = "csv")
# drive_download("~/tutorials/2019_03_variant_prioritization/ashkenazim.small_variants.csv")
# assignment <- Alt - 
id <- "15bYaktk50qD6qeL-LPgm5NSlk4-XFfYq"
# download by ID - sometimes you have to do it several times if it failed
drive_download(as_id(id), overwrite = T)
# download by URL
drive_download(as_id("https://drive.google.com/open?id=15bYaktk50qD6qeL-LPgm5NSlk4-XFfYq"), overwrite = T)

# some values are None in Quality and it reads as character string rather than numerical value
variants <- read_csv("ashkenazim.small_variants.csv", na = c("None"))
spec(variants) 
# remove low quality variants supported by 2 variant callers
variants <- filter(variants, Number_of_callers > 2)

# do it in excel first            
variants_rare_missense <- filter(variants, Variation == "missense_variant", Maf_all < 0.01)
variants_potentially_deleterious <- filter(variants_rare_missense, Polyphen_score > 0.9)

# explore variation
# pipe = Ctrl + Shift + M %>% 
variation <- variants %>% count(Variation)
# equivalent expression
# count(variants, Variation)

ggplot(data = variants) + geom_bar(mapping = aes(x = Variation)) + coord_flip()

variants_recessive <- filter(variants, Zygosity.Ashkenazim_HG002 == "Hom")

variants_recessive <- filter(variants, Zygosity.Ashkenazim_HG002 == "Hom", 
                                       Zygosity.Ashkenazim_HG003 == "Het",
                                       Zygosity.Ashkenazim_HG004 == "Het",
                                       Omim_gene_description != "0")
write_excel_csv(variants_recessive, "variants_recessive.csv")

# example of how to create and use a function
variants_recessive <- filter_variants(variants, "Hom", "Het", "Het") %>% 
    filter(Omim_gene_description != "0")

# X linked
variants_x_linked <- filter_variants(variants, "Hom","-", "Het") %>% 
    filter(grepl("X", Omim_inheritance))

# variants reported in clinvar 
clinvar <- variants %>% count(Clinvar)
# exersize : plot a bar chart

# pull only clinvar pathogenic variants
variants_clinvar_pathogenic <- filter(variants, Clinvar == "pathogenic")

# pull only variants with high CADD and polyphen scores
variants_computational <- filter(variants, Cadd_score > 30, Polyphen_score > 0.95, 
                        Zygosity.Ashkenazim_HG002 != "-", Omim_gene_description != "0")
# none of them is homozygous!
variants_computational <- filter(variants, Cadd_score > 30, Polyphen_score > 0.95, 
                        Zygosity.Ashkenazim_HG002 != "-", Omim_gene_description != "0",
                        Omim_inheritance = "AD")

# de novo
variants_de_novo <- filter_variants(variants, "Het", "-", "-")
variants_de_novo <- filter_variants(variants, "Het", "-", "-") %>% filter(Quality > 100)

# DNMT is a candate?
# save de novos
write_excel_csv(variants_de_novo, "variants_denovo.csv", append = F)

# Exercise: are there some more de novo variants?
variants_de_novo <- filter_variants(variants, "Hom", "Het", "-") %>% filter(Quality > 100)
variants_de_novo <- filter_variants(variants, "Hom", "-", "Het") %>% filter(Quality > 100)
write_excel_csv(variants_de_novo, "variants_denovo.csv", append = T)

# review the file - why so many de novo variants on chrX?
# is proband a male? - search for Ashkenazim trio
# (son - father - mother)
variants_chrx <- filter(variants, grepl("chrX", Position), Number_of_callers > 2)
# has a paralogue in chrY pseudoatosomal region
# https://www.ncbi.nlm.nih.gov/gene/414
variants_chrx <- filter(variants_chrx, Gene != "ARSD")
# what is wrong with ESX1 gene?
# low alt depth - false positive
# Het in ASMTL? - pseudoatosomal region (<2.7M)
# https://www.ncbi.nlm.nih.gov/grc/human
variants_chrx <- filter(variants_chrx, Gene != "ESX1")
# all these variants are not de novos!

# select a subset of columns
colnames(variants)
variants_de_novo_report <- select(variants_de_novo, c("Position", "Ref", "Alt", "Gene", "Variation",
                                                      "Gnomad_maf", "rsIDs"))
# exercise : select columns you'd want to report

# TCAG small variant report for WGS for NA12878
# https://drive.google.com/open?id=1RATcIzqoptbyfSoJHTxAUoHT23wCovcz
drive_download(as_id("https://drive.google.com/open?id=1RATcIzqoptbyfSoJHTxAUoHT23wCovcz"), overwrite = T)
wgs_small_variants_tcag <- read_tsv("NA12878.hc.vqsr.vcf_ext.tsv.annovar.out_SUBSET_rev26.5.chr22.tsv")

# pull variants matching HGMD database 
variants_hgmd <- filter(wgs_small_variants_tcag, grepl("DM",HGMD_tag))

# exercise: choose your strategy on how to extract pathogenic variants, filter + select with pipe %>% 

# structural variant report - for muscular genes
drive_download(as_id("https://drive.google.com/open?id=1O2RtjYbP2kyUxd-VAYhb3WPhgilJOmhj"), overwrite = T)
variants_sv <- read_csv("NA12878.sv.csv")

# sort by lengths decreasing
variants_sv <- arrange(variants_sv, SVLEN)
drive_download(as_id("https://drive.google.com/open?id=1WYPGlOMPvxITzKJg6YlWGajNfu92yyIP"), overwrite = T)
drive_download(as_id("https://drive.google.com/open?id=1WED_JAZJvRFKB3h7MqXSSfjqVT_1Ficr"), overwrite = T)
