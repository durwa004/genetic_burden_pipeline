library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(biomaRt)
library(gt)
library(biomartr)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(org.At.tair.db)
library(GOSim)
library(ggmosaic)
library(ggrepel)
library(webshot2)

setwd("/Users/durwa004/Documents/Research/GB_project")

high_int <- read.table("SnpEff.VEP.intersect.txt", header = T)
mean(high_int$AF)
range(high_int$AF)
length(high_int$Intersect)
table(high_int$Intersect)

#Get no. SNPs vs indels
snp_indels <- high_int %>% mutate(
  indel = nchar(alt) >1
)

table(snp_indels$indel)

#Get non GB AF
non <- read.table("non_intersect.txt", header=F)
mean(non$V1)
range(non$V1)

#Check for differences
t.test(high_int$AF, non$V1)

###Get count of genotypes for each individual
ind <- read.table("SnpEff.VEP.intersect.individual.txt", header = T)
ind$CHROM <- NULL
ind$POS <- NULL
ind$REF <- NULL
ind$ALT <- NULL

ind_ind <- ind %>% summarise_all( ~ list(table(.)))
#Based on output from python Get_genetic_burden_by_individual.py
het <- colSums(ind == "0/1") + colSums(ind == "0|1") + 
  colSums(ind == "0/2") + colSums(ind == "0|2") + 
  colSums(ind == "1/2") + colSums(ind == "1|2") + 
  colSums(ind == "0/3") + colSums(ind == "0|3") + 
  colSums(ind == "1/0") + colSums(ind == "1|0") + 
  colSums(ind == "1/3") + colSums(ind == "1|3") + 
  colSums(ind == "0/4") + colSums(ind == "0|4")
hom <- colSums(ind == "1/1") + colSums(ind == "1|1") +
  colSums(ind == "2/2") + colSums(ind == "2|2") + 
  colSums(ind == "3/3") + colSums(ind == "3|3")
mean(hom)
range(hom)
total_variant <- het + hom
mean(total_variant)
range(total_variant)

#Get number of genes
#Need to convert gene names so they're interpretable
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "ecaballus_gene_ensembl",
                      mirror = "useast")
searchAttributes(mart = ensembl, pattern = "ensembl_gene_id")

SE.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                values = high_int$SE.Gene_Name, 
                mart = ensembl)
VEP.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                     values = high_int$VEP.Gene, 
                     mart = ensembl)

SE.HGNC <- SE.gene.ID %>%
  subset(nchar(hgnc_symbol) >1)

#Merge with high_int
colnames(SE.gene.ID) <- c("hgnc_symbol", "SE.Gene_Name")
high_int_m <- merge(high_int, SE.gene.ID, by = "SE.Gene_Name")

gene_count <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  count(SE.Gene_Name)

length(gene_count$SE.Gene_Name)
mean(gene_count$n)
range(gene_count$n)
table(gene_count$n)

hgnc_count <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  count(hgnc_symbol)

hgnc_count <- hgnc_count[-1,]
length(hgnc_count$hgnc_symbol)
mean(hgnc_count$n)
range(hgnc_count$n)
table(hgnc_count$n > 5)

#Get AF in each gene
gene_count2 <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  summarize(AF = mean(AF))

mean(gene_count2$AF)
range(gene_count2$AF)

hgnc_count2 <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  summarize(AF = mean(AF))

mean(hgnc_count2$AF)
range(hgnc_count2$AF)

#Get genes containing 5 or more variants
high_hgnc <- hgnc_count %>% 
  filter(n >5)
length(high_hgnc$n)

#Export for DAVID analysis
hgnc_table <- as.data.frame(high_int_m$SE.Gene_Name)
hgnc_table$go_id <- high_int_m$hgnc_symbol

write.table(unique(hgnc_table$SE.Gene_Name), file = "hgnc_background.txt", quote = FALSE, sep = "\t")
write.table(high_hgnc$hgnc_symbol, file = "high_hgnc.txt", quote = FALSE, sep = "\t")

#Export for paper
high_hgnc <- as.data.frame(high_hgnc)
gt_high <- gt(high_hgnc)
gt_high <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 2. Count of HGNC gene IDs containing >5 genetic burden variants**")
  ) %>%
  cols_label(hgnc_symbol = md("**Gene ID**"),
             n = md("**Number**")) %>%
  cols_width(
    hgnc_symbol ~ px(150),
    n ~ px(100)
  )
gt_high %>% gtsave("supp_table2.pdf")

#Get genes with AF >50%
high_AF_hgnc <- hgnc_count2 %>%
  filter(AF >0.05)
length(high_AF_hgnc$AF)
write.table(high_AF_hgnc$hgnc_symbol, file = "high_AF_hgnc.txt", quote = FALSE, sep = "\t")

#Export for DAVID
high_hgnc$hgnc_symbol
high_hgnc <- merge(high_hgnc, SE.gene.ID, by = "hgnc_symbol")
high_AF_hgnc <- merge(high_AF_hgnc, SE.gene.ID, by = "hgnc_symbol")

#Export for paper
high_AF_hgnc <- as.data.frame(high_AF_hgnc)
gt_high_AF <- gt(high_AF_hgnc)
gt_high_AF <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 3. Count of HGNC gene IDs containing a genetic burden variants at an allele frequency >5%**")
  ) %>%
  cols_label(hgnc_symbol = md("**Gene ID**"),
             n = md("**Number**")) %>%
  cols_width(
    hgnc_symbol ~ px(150),
    n ~ px(100)
  )
gt_high_AF %>% gtsave("supp_table3.pdf")

##Frequency of loss of function variants
table(high_int$SE.LOF)
table(high_int$VEP.Consequence)

#Take first consequence when >1 consequence listed)
high_int_tidy <- high_int %>%
  mutate(across('VEP.Consequence', str_replace, "&inframe_deletion", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&splice_region_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&3_prime_UTR_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&splice_region_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&protein_altering_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&inframe_insertion", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&inframe_deletion", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&frameshift_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&start_retained_variant", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&5_prime_UTR_variant", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&stop_lost", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&splice_acceptor_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&intron_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&coding_sequence_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&non_coding_transcript_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&missense_variant", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&non_coding_transcript_exon_variant", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&splice_donor_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "&stop_retained_variant", "")) %>% 
  mutate(across('VEP.Consequence', str_replace, "&start_lost", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&inframe_deletion", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&splice_region_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&3_prime_UTR_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&splice_region_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&protein_altering_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&inframe_insertion", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&inframe_deletion", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&frameshift_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&start_retained_variant", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&5_prime_UTR_variant", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&stop_lost", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&splice_acceptor_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&intron_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&coding_sequence_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&non_coding_transcript_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&missense_variant", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&non_coding_transcript_exon_variant", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&splice_donor_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&stop_retained_variant", "")) %>% 
  mutate(across('SE.Annotation', str_replace, "&start_lost", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&disruptive_inframe_insertion", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&disruptive_inframe_deletion", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&conservative_inframe_insertion", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&conservative_inframe_deletion", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&3_prime_UTR_truncation", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&exon_loss_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "&stop_gained", "")) %>%
  mutate(across('SE.Annotation', str_replace, "conservative_inframe_deletion", "inframe deletion")) %>%
  mutate(across('SE.Annotation', str_replace, "conservative_inframe_insertion", "inframe_insertion")) %>%
  mutate(across('SE.Annotation', str_replace, "disruptive_inframe_deletion", "inframe_deletion")) %>%
  mutate(across('SE.Annotation', str_replace, "_variant", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "_variant", "")) %>%
  mutate(across('SE.Annotation', str_replace, "_", " ")) %>%
  mutate(across('VEP.Consequence', str_replace, "_", " ")) %>%
  mutate(across('SE.Annotation', str_replace, "_", " ")) %>%
  mutate(across('VEP.Consequence', str_replace, "_", " ")) %>%
  mutate(across('SE.Annotation', str_replace, "bidirectional ", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "bidirectional", ""))

high_int_tidy <- high_int_tidy %>%
  mutate(consequence_agreement = if_else(SE.Annotation == VEP.Consequence, "yes", "no"))

table(high_int_tidy$consequence_agreement)
      
#x <- ggplot(data = high_int_tidy) + theme_bw() + 
#  geom_mosaic(aes(x = product(SE.Annotation, VEP.Consequence), 
#                  fill = consequence_agreement)) +   
#  labs(y="Ensembl VEP consequence", x="SnpEff annotation") + 
#  guides(x = guide_axis(n.dodge = 15)) + 
#  theme(panel.grid = element_blank(), panel.border = element_blank(), 
#        axis.line.x = element_line(), axis.line.y = element_line(), 
#        axis.text = element_text(size=18), 
#        axis.title = element_text(size=16,face="bold"),
#        legend.position="none")
#save_plot("SE_VEP_consequence.tiff", x, base_height = 18, base_width = 18)
#x <- ggplot(data = high_int_tidy) + theme_bw() + 
#  geom_mosaic(aes(x = product(SE.Annotation), 
#                  fill = VEP.Consequence)) +   
#  geom_mosaic_text(aes(x = product(SE.Annotation, VEP.Consequence)), 
#                   stat = "mosaic", position = "identity", offset = 0.01,
#                   show.legend = NA, repel = TRUE, as.lab = TRUE) + 
#  labs(y="Ensembl VEP consequence", x="SnpEff annotation", 
#       title = "Figure 2. Variant consequence determined by Ensembl VEP and SnpEff",
#       subtitle = "Mosaic plot of variant consequence overlap between Ensembl VEP and SnpEff. Color based on\nEnsembl VEP consequence") + 
#  theme(panel.grid = element_blank(), panel.border = element_blank(), 
#        axis.text.x = element_blank(),
#        axis.text.y = element_blank(),
#        axis.text = element_text(size=8), axis.title = element_text(size=10,face="bold"),
#        legend.title = element_blank(),
#        plot.title = element_text(size=10, face = "bold", hjust = 0),
#        plot.subtitle = element_text(size=8, hjust = 0),
#        plot.title.position = "plot")
#save_plot("figure2.tiff", x, base_height = 8, base_width = 6)

##Get table of consequences for SnpEff and VEP
#Need to get a list of all the consequences of VEP and SnpEff 
#and then get counts for each one

#high_int_tidy$consequence_combined <- paste(high_int_tidy$SE.Annotation,
                                            #high_int_tidy$VEP.Consequence, sep = "_")
  
consequence <- c(high_int_tidy$SE.Annotation, high_int_tidy$VEP.Consequence)
consequence <- unique(consequence)

SE_c <- as.table(colSums(sapply(consequence, grepl, high_int_tidy$SE.Annotation)))
VEP_c <- as.table(colSums(sapply(consequence, grepl, high_int_tidy$VEP.Consequence)))

SE_VEP_tbl <- as.data.frame(cbind(consequence, SE_c, VEP_c))
SE_VEP_tbl$SE_c <- as.numeric(SE_VEP_tbl$SE_c)
SE_VEP_tbl$VEP_c <- as.numeric(SE_VEP_tbl$VEP_c)

gt_SE_VEP <- gt(SE_VEP_tbl)
gt_SE_VEP <- 
  gt_SE_VEP %>%
  tab_header(
    title = md("**Table 2. Overlap between SnpEff and Ensembl-VEP\npredicted variant consequence**")
  )  %>%
  cols_label(consequence = md("**Consequence**"),
             SE_c = md("**SnpEff**"),
             VEP_c = md("**Ensembl VEP**")) %>%
  fmt_number(columns = 2:3, sep_mark = ",", use_seps = TRUE, decimals = 0) %>%
  cols_width(
    consequence ~ px(150),
    SE_c ~ px(100),
    VEP_c ~ px(100)
  )
gt_SE_VEP %>% gtsave("table2.pdf")


#Get LOF variants for SnpEff and Ensembl VEP 
##Need to figure it out manually for VEP
table(high_int_tidy$SE.LOF) #44581

high_int <- high_int %>%
  mutate(VEP.LOF = str_detect(VEP.Consequence, "frameshift|splice_acceptor_variant|
                              splice_donor_variant|splice_region|stop_gained"))
cons <- c("Not_LOF", "LOF")
high_int <- high_int %>%
  mutate(across("VEP.LOF", str_replace, "TRUE", "yes")) %>%
  mutate(across("VEP.LOF", str_replace, "FALSE", "no"))

high_int$LOF <- paste(high_int$SE.LOF, high_int$VEP.LOF, sep = ":")
table(high_int$LOF)

#LOF per individual
ind <- read.table("SnpEff.VEP.intersect.individual.txt", header = T)
ind$chrom_pos <- paste(ind$CHROM, ind$POS, sep = ":")
high_int$chrom_pos <- paste(high_int$chrom, high_int$pos, sep = ":")

LOF <- high_int$chrom_pos[high_int$LOF == "yes:yes"]
mean(high_int$AF[high_int$LOF == "yes:yes"])
range(high_int$AF[high_int$LOF == "yes:yes"])

LOF_ind <- ind %>% 
  filter(ind$chrom_pos %in% LOF)
LOF_ind$CHROM <- NULL
LOF_ind$POS <- NULL
LOF_ind$REF <- NULL
LOF_ind$ALT <- NULL
LOF_ind$chrom_pos <- NULL

LOF_ind_ind <- LOF_ind %>% summarise_all( ~ list(table(.)))
#Based on output from python Get_genetic_burden_by_individual.py
het <- colSums(LOF_ind == "0/1") + colSums(LOF_ind == "0|1") + 
  colSums(LOF_ind == "0/2") + colSums(LOF_ind == "0|2") + 
  colSums(LOF_ind == "1/2") + colSums(LOF_ind == "1|2") + 
  colSums(LOF_ind == "0/3") + colSums(LOF_ind == "0|3") + 
  colSums(LOF_ind == "1/0") + colSums(LOF_ind == "1|0") + 
  colSums(LOF_ind == "1/3") + colSums(LOF_ind == "1|3") + 
  colSums(LOF_ind == "0/4") + colSums(LOF_ind == "0|4")
hom <- colSums(LOF_ind == "1/1") + colSums(LOF_ind == "1|1") +
  colSums(LOF_ind == "2/2") + colSums(LOF_ind == "2|2") + 
  colSums(LOF_ind == "3/3") + colSums(LOF_ind == "3|3")
mean(hom)
range(hom)
total_variant <- het + hom
mean(total_variant)
range(total_variant)


#Get LOF genes
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "ecaballus_gene_ensembl",
                      mirror = "useast")
searchAttributes(mart = ensembl, pattern = "ensembl_gene_id")

SE.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                     values = high_int$SE.Gene_Name[high_int$LOF == "yes:yes"], 
                     mart = ensembl)
VEP.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                      values = high_int$VEP.Gene[high_int$LOF == "yes:yes"], 
                      mart = ensembl)

SE.HGNC <- SE.gene.ID %>%
  subset(nchar(hgnc_symbol) >1)

#Merge with high_int
colnames(SE.gene.ID) <- c("hgnc_symbol", "SE.Gene_Name")
high_int_m <- merge(high_int[high_int$LOF == "yes:yes",], SE.gene.ID, by = "SE.Gene_Name")

gene_count <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  count(SE.Gene_Name)

length(unique(gene_count$SE.Gene_Name))
mean(gene_count$n)
range(gene_count$n)
table(gene_count$n)

hgnc_count <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  count(hgnc_symbol)

#Remove ""
hgnc_count <- hgnc_count[-1,]
length(hgnc_count$hgnc_symbol)
mean(hgnc_count$n)
range(hgnc_count$n)
table(hgnc_count$n)
table(hgnc_count$n > 5)

#Get AF in each gene
gene_count2 <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  summarize(AF = mean(AF))

mean(gene_count2$AF)
range(gene_count2$AF)

hgnc_count2 <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  summarize(AF = mean(AF))

mean(hgnc_count2$AF)
range(hgnc_count2$AF)

#Get genes containing 5 or more variants
high_hgnc <- hgnc_count %>% 
  filter(n >5)
length(high_hgnc$n)

#Export for DAVID analysis
hgnc_table <- as.data.frame(high_int_m$SE.Gene_Name)
hgnc_table$go_id <- high_int_m$hgnc_symbol
colnames(hgnc_table) <- c("SE.Gene_Name", "go_id")
write.table(unique(hgnc_table$SE.Gene_Name), file = "LOF_hgnc_background.txt", quote = FALSE, sep = "\t")
write.table(high_hgnc$hgnc_symbol, file = "LOF_high_hgnc.txt", quote = FALSE, sep = "\t")

#Export for paper
high_hgnc <- as.data.frame(high_hgnc)
gt_high <- gt(high_hgnc)
gt_high <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 4. Count of HGNC gene IDs containing >5 LOF variants**")
  ) %>%
  cols_label(hgnc_symbol = md("**Gene ID**"),
             n = md("**Number**")) %>%
  cols_width(
    hgnc_symbol ~ px(150),
    n ~ px(100)
  )
gt_high %>% gtsave("supp_table4.pdf")

#Get genes with AF >50%
high_AF_hgnc <- hgnc_count2 %>%
  filter(AF >0.05)
length(high_AF_hgnc$AF)
write.table(high_AF_hgnc$hgnc_symbol, file = "LOF_high_AF_hgnc.txt", quote = FALSE, sep = "\t")

#Export for DAVID
high_hgnc$hgnc_symbol
high_hgnc <- merge(high_hgnc, SE.gene.ID, by = "hgnc_symbol")
high_AF_hgnc <- merge(high_AF_hgnc, SE.gene.ID, by = "hgnc_symbol")

#Export for paper
high_AF_hgnc <- as.data.frame(high_AF_hgnc)
gt_high_AF <- gt(high_AF_hgnc)
gt_high_AF <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 5. Count of HGNC gene IDs containing LOF variants at an allele frequency >5%**")
  ) %>%
  cols_label(hgnc_symbol = md("**Gene ID**"),
             n = md("**Number**")) %>%
  cols_width(
    hgnc_symbol ~ px(150),
    n ~ px(100)
  )
gt_high_AF %>% gtsave("supp_table5.pdf")


#Look at hand annotated variants to see overlap - compare to high_int_tidy - there is no overlap!
hand <- read.table("hand_validated_tidy.txt", header=TRUE)
hand <- hand %>%
  mutate(across('chrom', str_replace, "NC_009144.3", "chr1")) %>%
  mutate(across('chrom', str_replace, "NC_009145.3", "chr2")) %>%
  mutate(across('chrom', str_replace, "NC_009146.3", "chr3")) %>%
  mutate(across('chrom', str_replace, "NC_009147.3", "chr4")) %>%
  mutate(across('chrom', str_replace, "NC_009148.3", "chr5")) %>%
  mutate(across('chrom', str_replace, "NC_009149.3", "chr6")) %>%
  mutate(across('chrom', str_replace, "NC_009150.3", "chr7")) %>%
  mutate(across('chrom', str_replace, "NC_009151.3", "chr8")) %>%
  mutate(across('chrom', str_replace, "NC_009152.3", "chr9")) %>%
  mutate(across('chrom', str_replace, "NC_009153.3", "chr10")) %>%
  mutate(across('chrom', str_replace, "NC_009154.3", "chr11")) %>%
  mutate(across('chrom', str_replace, "NC_009155.3", "chr12")) %>%
  mutate(across('chrom', str_replace, "NC_009156.3", "chr13")) %>%
  mutate(across('chrom', str_replace, "NC_009157.3", "chr14")) %>%
  mutate(across('chrom', str_replace, "NC_009158.3", "chr15")) %>%
  mutate(across('chrom', str_replace, "NC_009159.3", "chr16")) %>%
  mutate(across('chrom', str_replace, "NC_009160.3", "chr17")) %>%
  mutate(across('chrom', str_replace, "NC_009161.3", "chr18")) %>%
  mutate(across('chrom', str_replace, "NC_009162.3", "chr19")) %>%
  mutate(across('chrom', str_replace, "NC_009163.3", "chr20")) %>%
  mutate(across('chrom', str_replace, "NC_009164.3", "chr21")) %>%
  mutate(across('chrom', str_replace, "NC_009165.3", "chr22")) %>%
  mutate(across('chrom', str_replace, "NC_009166.3", "chr23")) %>%
  mutate(across('chrom', str_replace, "NC_009167.3", "chr24")) %>%
  mutate(across('chrom', str_replace, "NC_009168.3", "chr25")) %>%
  mutate(across('chrom', str_replace, "NC_009169.3", "chr26")) %>%
  mutate(across('chrom', str_replace, "NC_009170.3", "chr27")) %>%
  mutate(across('chrom', str_replace, "NC_009171.3", "chr28")) %>%
  mutate(across('chrom', str_replace, "NC_009172.3", "chr29")) %>%
  mutate(across('chrom', str_replace, "NC_009173.3", "chr30")) %>%
  mutate(across('chrom', str_replace, "NC_009174.3", "chr31")) %>%
  mutate(across('chrom', str_replace, "NC_009175.3", "chrX"))

hand$chrom_pos <- paste(hand$chrom, hand$pos, sep = ":")

high_int_tidy$chrom_pos <- paste(high_int_tidy$chrom, high_int_tidy$pos, sep = ":")
high_int_hand <- high_int_tidy %>%
  mutate(hand = if_else(high_int_tidy$chrom_pos == hand$chrom_pos, "yes", "no"))

high_int_hand <- high_int_tidy$chrom_pos[hand$chrom_pos]
high_int_hand <- as.factor(high_int_hand)

###No gene enrichment using enricher
#test <- biomartr::organismBM(organism = "Equus caballus")
ecaballus_attributes <- 
  biomartr::organismAttributes("Equus caballus") %>% 
  filter(dataset == "ecaballus_gene_ensembl")

attributes_to_retrieve = c("hgnc_id", "entrezgene_id", "go_id", "name_1006", "namespace_1003")

result_BM <- biomartr::biomart( genes      = high_AF_hgnc$SE.Gene_Name,                  # genes were retrieved using biomartr::getGenome()
                                mart       = "ENSEMBL_MART_ENSEMBL",                     # marts were selected with biomartr::getMarts()
                                dataset    = "ecaballus_gene_ensembl",               # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
                                filters =   "ensembl_gene_id" )# query key
head(result_BM)  
all_EC_annotated <- biomartr::biomart(genes = high_int$SE.Gene_Name,
                                      mart       = "ENSEMBL_MART_ENSEMBL", 
                                      dataset    = "ecaballus_gene_ensembl",  
                                      attributes = attributes_to_retrieve,
                                      filters =  "ensembl_gene_id" )

EC_t2g <- as.data.frame(result_BM$go_id)
EC_t2g$go_di <- result_BM$name_1006

ora_analysis_high <- enricher(gene = result_BM$name_1006, 
                            universe = all_EC_annotated$name_1006, 
                            pAdjustMethod = "bonferroni",
                            minGSSize = 2,
                            maxGSSize = 500,
                            TERM2GENE = EC_t2g,
                            pvalueCutoff = 0.05)

dotplot(ora_analysis_high)

ora_analysis_high <- as.data.frame(ora_analysis_high)
ora_analysis_bp <- pairwise_termsim(ora_analysis_high, method = "JC")
emapplot(ora_analysis_high, color = "qvalue")



