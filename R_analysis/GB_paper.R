library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)
library(devtools)
library(ggpubr)
library(scales)
library(emmeans)
library(dvmisc)
library(tidyr)
library(RColorBrewer)
library(tidyverse)
library(gt)
library(biomaRt)
library(gtsummary)
library(webshot2)
library(pals)
library(ggpattern)

setwd("/Users/durwa004/Documents/Research/GB_project/")

####Equine genetic variation catalog summary####

##Table 1 - horse breeds and DOC
breeds <- read.table("horse_genomes_breeds_all.txt", header=T)
colnames(breeds) <- c("ID", "Breed")
breeds <- breeds[!(breeds$ID %in% c("M989","M6468")),]

#Get DOC
bcftools <- read.table("subset_number_of_variants.txt", header=T)
bcftools$nVariants <- bcftools$nSNPs + bcftools$indels

bcf_breed <- merge(breeds, bcftools, by = "ID")

bcf_breed <- bcf_breed %>%
  mutate(across('Breed', str_replace_all, "_", " ")) %>%
  mutate(across('Breed', str_replace_all, "_x", " x")) %>%
  mutate(across('Breed', str_replace, "Uknown", "Unknown")) %>%
  mutate(across('Breed', str_replace, "Trakenher", "Trakehner")) %>%
  mutate(across('Breed', str_replace, "Halflinger", "Haflinger")) %>%
  mutate(across('Breed', str_replace, "Hanovarian", "Hanoverian")) %>%
  mutate(across('Breed', str_replace, "ColdBlood", "Coldblood")) %>%
  mutate(across('Breed', str_replace, "Other horse", "Unknown")) %>%
  mutate(across('Breed', str_replace, "Pony", "Unknown")) %>%
  mutate(across('Breed', str_replace, "QH", "Quarter Horse")) %>%
  mutate(across('Breed', str_replace, "STB", "Standardbred")) %>%
  mutate(across('Breed', str_replace, "TB", "Thoroughbred")) %>%
  mutate(across('Breed', str_replace, "WarmBlood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "TWH", "Tennessee Walking Horse")) %>%
  mutate(across('Breed', str_replace, "duelmener", "Dulmener"))

mean(bcf_breed$DOC)
range(bcf_breed$DOC)
mean(bcf_breed$nVariants)
mean(bcf_breed$nSNPs)
mean(bcf_breed$indels)
mean(bcf_breed$nSingletons)

breed_tbl <- bcf_breed %>%
  group_by(Breed) %>%
  summarize(count=n(), m=mean(DOC))
breed_tbl$m <- lapply(breed_tbl$m, round, 1)

colnames(breed_tbl) <- c("Breed", "Number", "DOC")
breed_tbl <- breed_tbl %>% arrange(Breed)

gt_breed <- gt(breed_tbl)
gt_breed <- 
  gt_breed %>%
  tab_header(
    title = md("**Supplementary table 1.** Breed information for the 605 horses")
  ) %>%
  cols_label(Breed = md("**Breed**"),
             Number = md("**Number**"),
             DOC = md("**Average DOC**")) %>%
  tab_footnote(
    footnote = "Breeds with 17 or more horses.",
    locations = cells_body(
      columns = Number,
      rows = Number >=17
    )) %>% 
  cols_width(
    Breed ~ px(150),
    Number ~ px(100),
    DOC ~ px(100)
  ) %>%   tab_options(
    table.font.size = 8,
    table.width = px(200)) 
  )

gt_breed %>% gtsave("Figures_tables/supp_table1.pdf", expand = 0)


###Number of variants by individual
bcftools_br <- bcf_breed %>%
  mutate(across('Breed', str_replace_all, "_", " ")) %>%
  mutate(across('Breed', str_replace, "Uknown", "Other")) %>%
  mutate(across('Breed', str_replace, "Yakut", "Other")) %>%
  mutate(across('Breed', str_replace, "Westphalian", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Trakenher", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Unknown", "Other")) %>%
  mutate(across('Breed', str_replace, "UK Warmblood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Trakehner", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Swiss Warmblood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Sports Horse", "Other")) %>%
  mutate(across('Breed', str_replace, "Sorraia", "Other")) %>%
  mutate(across('Breed', str_replace, "Saxon-Thuringian Heavy Warmblood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Saddle Trotter", "Other")) %>%
  mutate(across('Breed', str_replace, "Percheron", "Other")) %>%
  mutate(across('Breed', str_replace, "Paint", "Quarter Horse")) %>%
  mutate(across('Breed', str_replace, "Oldenberg", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Norwegian Fjord", "Other")) %>%
  mutate(across('Breed', str_replace, "Native Mongolian Chakouyi Horse", "Other")) %>%
  mutate(across('Breed', str_replace, "Mongolian", "Other")) %>%
  mutate(across('Breed', str_replace, "Missouri Fox Trotter x", "Other")) %>%
  mutate(across('Breed', str_replace, "Miniature Horse", "Other")) %>%
  mutate(across('Breed', str_replace, "Mangalarga Marchador Horse", "Other")) %>%
  mutate(across('Breed', str_replace, "Lipizzaner", "Other")) %>%
  mutate(across('Breed', str_replace, "KWPN", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Italian trotter", "Other")) %>%
  mutate(across('Breed', str_replace, "Holsteiner", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Hanoverian", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "German WB", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Friesian", "Other")) %>%
  mutate(across('Breed', str_replace, "French Trotter", "Other")) %>%
  mutate(across('Breed', str_replace, "Draft x", "Other")) %>%
  mutate(across('Breed', str_replace, "Dulmener", "Other")) %>%
  mutate(across('Breed', str_replace, "Curly Trotter", "Other")) %>%
  mutate(across('Breed', str_replace, "Connemara", "Other")) %>%
  mutate(across('Breed', str_replace, "Coldblood", "Other")) %>%
  mutate(across('Breed', str_replace, "British Warmblood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Aegidienberger Paso x Icelandic", "Other")) %>%
  mutate(across('Breed', str_replace, "Halflinger", "Other")) %>%
  mutate(across('Breed', str_replace, "Haflinger", "Other")) %>%
  mutate(across('Breed', str_replace, "Hanovarian", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "ColdBlood", "Other")) %>%
  mutate(across('Breed', str_replace, "Other horse", "Other")) %>%
  mutate(across('Breed', str_replace, "Pony", "Other")) %>%
  mutate(across('Breed', str_replace, "QH", "Quarter Horse")) %>%
  mutate(across('Breed', str_replace, "STB", "Standardbred")) %>%
  mutate(across('Breed', str_replace, "TB", "Thoroughbred")) %>%
  mutate(across('Breed', str_replace, "WarmBlood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "TWH", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "duelmener", "Other")) %>%
  mutate(across('Breed', str_replace, "UK Warmblood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "Tennessee Walking Horse", "Other")) %>%
  mutate(across('Breed', str_replace, "WP", "Welsh Pony")) %>%
  mutate(across('Breed', str_replace, "Jeju pony", "Other")) %>%
  mutate(across('Breed', str_replace, "POA", "Other")) %>%
  mutate(across('Breed', str_replace, "Quarter Horse x", "Other")) %>%
  mutate(across('Breed', str_replace, "Jeju pony", "Other"))
  mutate(across('Breed', str_replace, "Thoroughbred x", "Other"))

bcftools_br$br[bcftools_br$Breed == "Arabian"] <- "Arab"
bcftools_br$br[bcftools_br$Breed == "Belgian"] <- "Belg"
bcftools_br$br[bcftools_br$Breed == "Clydesdale"] <- "Clyd"
bcftools_br$br[bcftools_br$Breed == "Franchese Montagne"] <- "FM"
bcftools_br$br[bcftools_br$Breed == "Icelandic"] <- "Ice"
bcftools_br$br[bcftools_br$Breed == "Morgan"] <- "Morg"
bcftools_br$br[bcftools_br$Breed == "Other"] <- "Oth"
bcftools_br$br[bcftools_br$Breed == "Quarter Horse"] <- "QH"
bcftools_br$br[bcftools_br$Breed == "Shetland"] <- "Shet"
bcftools_br$br[bcftools_br$Breed == "Standardbred"] <- "STB"
bcftools_br$br[bcftools_br$Breed == "Thoroughbred"] <- "TB"
bcftools_br$br[bcftools_br$Breed == "Warmblood"] <- "WB"
bcftools_br$br[bcftools_br$Breed == "Welsh Pony"] <- "WP"

color_mapping <- c("Arab" = "grey", "Belg" = "blue", "Clyd" = "cyan",
                   "FM" = "red", "Ice" = "purple", "Morg" = "orchid1",
                   "Oth" = "khaki", "QH" = "dodgerblue", "Shet" = "green4",
                   "STB" = "black", "TB" = "maroon", "WB" = "#FDBF6F",
                   "WP" = "yellow")


#Look for correlation between number of variants and DOC with line of best fit 
x = ggplot(bcftools_br, aes(x=DOC,y=nVariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Depth of coverage") + geom_point(aes(color=Breed)) + 
  scale_color_manual(values = as.vector(polychrome(13))) +
  scale_x_continuous(breaks = seq(0,50,5))+
  scale_y_continuous(breaks = seq(0,8000000,1000000), labels=comma, limits = c(0,8000000)) + geom_smooth() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), 
        axis.text = element_text(size=8), axis.title = element_text(size=10,face="bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(1, "mm"),
        plot.title = element_text(size=10, face = "bold", hjust = 0),
        plot.subtitle = element_text(size=8, hjust = 0),
        plot.title.position = "plot")

save_plot("Figures_tables/figure_1.tiff", x, base_height = 4, base_width = 6)


#Show interaction between nvariants and DOC
fit1 <- (lm(nVariants ~ Breed, data=bcftools_br))
fit2 <- (lm(nVariants ~ Breed + DOC, data=bcftools_br))
anova(fit1,fit2)

variants_DOC <- (lm(nVariants ~ br + DOC,data=bcftools_br))
summary(variants_DOC)

DOC_br <- (lm(DOC ~ br,data=bcftools_br))
summary(DOC_br)

#https://github.com/rvlenth/emmeans/blob/master/R/plot.emm.R
#Get EMMEANs
br_emm <- emmeans(DOC_br, specs = "br", weights = "proportional", 
                  type = "response")
br_emm2 <- as.data.frame(br_emm)
br_emm2 <- br_emm2 %>%
  arrange(emmean)

DOC_emm <- emmeans(variants_DOC, specs = "br", weights = "proportional", 
                  type = "response")
DOC_emm2 <- as.data.frame(DOC_emm)
DOC_emm2 <- DOC_emm2 %>%
  arrange(emmean)


x1 <- plot(br_emm2) + geom_boxplot(col = color_mapping[br_emm2$br]) + theme_bw() + 
  xlab("EMMEAN of depth of coverage") + ylab("Breed") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text = element_text(size=), axis.title = element_text(size=12,face="bold"))

x2 <- plot(DOC_emm2) + geom_boxplot(col = color_mapping[DOC_emm2$br]) + 
  theme_bw() + xlab("EMMEAN of Number of variants") + ylab("Breed") +
  scale_x_continuous(labels=comma) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), 
        axis.text = element_text(size=8), axis.title = element_text(size=12,face="bold"))

c <- ggarrange(x, ggarrange(x1, x2, ncol = 2, labels = c("b", "c"), align = "h"),
               nrow = 2, labels = "a")
save_plot("Figures_tables/figure_1.pdf", c, base_height = 8, base_width = 8)


####Overlap between Ensembl-VEP and SnpEff impact####
##Table 1: Get count of overlap
int <- read.table("SnpEff.VEP.intersect.txt", header = T)
int$intersect <- paste(int$SE.Annotation_Impact, ":", int$VEP.IMPACT)
table(int$intersect)

int_high <- int %>%
  filter(intersect == "HIGH : HIGH")

int_high_mod <- int %>%
  filter(intersect == "HIGH : MODERATE") 

int_mod_high  <- int %>%
  filter(intersect == "MODERATE : HIGH")

int_all <- rbind(int_high,int_high_mod,int_mod_high)
table(int_all$intersect)
length(int_all$intersect)
length(which(int_all$alt == "A")) + length(which(int_all$alt == "G")) +
  length(which(int_all$alt == "C")) + length(which(int_all$alt == "T"))

#Keep only high:high, high:mod, or mod:high variants

int_tbl <- int %>%
  group_by(intersect) %>%
  summarize(count=n())
colnames(int_tbl) <- c("Impact", "Count")  


overlap_table <- table(int$SE.Annotation_Impact, int$VEP.IMPACT)
overlap_df <- as.data.frame(overlap_table)
overlap_df$Var1 <- factor(overlap_df$Var1, levels = c("HIGH", "MODERATE", "LOW"))  # Change the order as needed
overlap_df$Var2 <- factor(overlap_df$Var2, levels = c("HIGH", "MODERATE", "LOW"))  # Change the order as needed


# Create the heatmap-like table using ggplot2
heatmap_table <- ggplot(overlap_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::comma(Freq)), color = "black") +  # Add frequency values
  scale_fill_gradient(low = "white", high = "lightblue") +
  xlab("Ensembl-VEP impact") + ylab("SnpEff impact") + 
  labs(title = "Table 1. Overlap between SnpEff and Ensembl-VEP variant impact") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), 
        axis.text = element_text(size=8), axis.title = element_text(size=10,face="bold"),
        legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(size=10, face = "bold", hjust = 0),
        plot.subtitle = element_text(size=8, hjust = 0),
        plot.margin = margin(0.5, 0.5, 0.5, 05, "cm"),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),     # Remove strip background
        strip.text.x = element_blank(),         # Remove x-axis strip text
        strip.text.y = element_blank(),          # Remove y-axis strip text
        plot.title.position = "plot")+
  coord_fixed(ratio = 0.3) + 
  scale_x_discrete(expand = c(0, 0)) +  # Adjust x-axis expand
  scale_y_discrete(expand = c(0, 0)) 

# Display the heatmap-like table
print(heatmap_table)

ggsave("Figures_tables/table1.pdf", heatmap_table)


####Overlap between Ensembl-VEP and SnpEff impact####
##Table 1: Get count of overlap
int <- read.table("SnpEff_VEP_intersect.txt", header = F)
colnames(int) <- c("location", "SE_VEP")

int_tbl <- int %>%
  group_by(SE_VEP) %>%
  summarize(count=n())

int_tbl<- int_tbl %>%
  mutate(across('SE_VEP', str_replace, "_", ":"))
colnames(int_tbl) <- c("Impact", "Count")  


gt_int <- gt(int_tbl)
gt_int <- 
  gt_int %>%
  tab_header(
    title = md("**Table 1.** Overlap between SnpEff and Ensembl-VEP variant impact")
  )  %>%
  cols_label(Impact = md("**SnpEff:VEP impact**"),
             Count = md("**Count**")) %>%
  tab_style(
    style = cell_fill(color = "lightgray"),
    locations = cells_body(
      rows = Impact == "HIGH:HIGH")) %>%
  tab_style(
    style = cell_fill(color = "lightgray"),
    locations = cells_body(
      rows = Impact == "HIGH:MODERATE")) %>%
  tab_style(
    style = cell_fill(color = "lightgray"),
    locations = cells_body(
      rows = Impact == "MODERATE:HIGH")) %>%
  fmt_number(columns = Count, sep_mark = ",", use_seps = TRUE, decimals = 0)  %>% tab_options(
    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(200)) %>%
  cols_width(
    Impact ~ px(200),
    Count ~ px(100)
  )

gt_int %>% gtsave("Figures_tables/table1.png")

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
non <- read.table("non_intersect.txt", header=T)
summary(non$AF)

#Check for differences
t.test(int_all$AF, non$V1)
non <- read.table("non_intersect.txt", header=F)
mean(non$V1)
range(non$V1)

#Check for differences
t.test(high_int$AF, non$V1)


###Get count of genotypes for each individual
ind <- read.table("SnpEff.VEP.intersect.individual.txt", header = T)
#ind$chrom_pos <- paste(ind$CHROM, ind$POS, sep = ":")
#tmp <- merge(ind, high_int, by = "chrom_pos")

ind_details <- ind[,0:4]
ind$CHROM <- NULL
ind$POS <- NULL
ind$REF <- NULL
ind$ALT <- NULL

ind <- ind %>%
  mutate_all(., str_replace_all, "/", "") %>%
  mutate_all(., str_replace_all, "\\|", "") %>%
  mutate_all(., str_replace, "01", "het") %>%
  mutate_all(., str_replace, "02", "het") %>%
  mutate_all(., str_replace, "12", "het") %>%
  mutate_all(., str_replace, "12", "het") %>%
  mutate_all(., str_replace, "03", "het") %>%
  mutate_all(., str_replace, "03", "het") %>%
  mutate_all(., str_replace, "10", "het") %>%
  mutate_all(., str_replace, "10", "het") %>%
  mutate_all(., str_replace, "13", "het") %>%
  mutate_all(., str_replace, "13", "het") %>%
  mutate_all(., str_replace, "04", "het") %>%
  mutate_all(., str_replace, "04", "het") %>%
  mutate_all(., str_replace, "11", "hom") %>%
  mutate_all(., str_replace, "11", "hom") %>%
  mutate_all(., str_replace, "22", "hom") %>%
  mutate_all(., str_replace, "22", "hom") %>%
  mutate_all(., str_replace, "33", "hom") %>%
  mutate_all(., str_replace, "33", "hom") %>%
  mutate_all(., str_replace, "00", "hom_WT")  

#Get het/hom by individual
het <- colSums(ind == "het")
hom <- colSums(ind == "hom")
median(hom)
summary(hom)
total_variant <- het + hom
median(total_variant)
summary(total_variant)
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

colnames(SE.gene.ID) <- c("hgnc_symbol", "SE.Gene_Name", "N1", "N2")
high_int_m <- merge(int_all, SE.gene.ID, by = "SE.Gene_Name")


colnames(SE.gene.ID) <- c("hgnc_symbol", "SE.Gene_Name")
high_int_m <- merge(high_int, SE.gene.ID, by = "SE.Gene_Name")


gene_count <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  count(SE.Gene_Name)

length(unique(gene_count$SE.Gene_Name))

summary(gene_count$n)
table(gene_count$n)

gene_count2 <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  summarize(AF = mean(AF))
summary(gene_count2$AF)

hgnc_count <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  count(hgnc_symbol)

hgnc_count <- hgnc_count[-1,]
length(hgnc_count$hgnc_symbol)
summary(hgnc_count$n)
range(hgnc_count$n)
table(hgnc_count$n)
table(hgnc_count$n > 5)

hgnc_count2 <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  summarize(AF = mean(AF))

summary(hgnc_count2$AF)

mean(gene_count$n)
range(gene_count$n)
table(gene_count$n)


#hgnc_count <- high_int_m %>%
#  group_by(hgnc_symbol) %>%
#  count(hgnc_symbol)

#hgnc_count <- hgnc_count[-1,]
#length(hgnc_count$hgnc_symbol)
#mean(hgnc_count$n)
#range(hgnc_count$n)
#table(hgnc_count$n)
#table(hgnc_count$n > 5)

#hgnc_count2 <- high_int_m %>%
#  group_by(hgnc_symbol) %>%
#  summarize(AF = mean(AF))

#mean(hgnc_count2$AF)
#range(hgnc_count2$AF)


#Get genes containing 5 or more variants
#high_hgnc <- hgnc_count %>% 
#  filter(n >5)
#length(high_hgnc$n)

#Export for DAVID analysis
#hgnc_table$go_id <- high_int_m$hgnc_symbol

gene_count_h <- gene_count %>%
  filter(n > 5)
length(unique(gene_count_h$SE.Gene_Name))
write.table(unique(gene_count$SE.Gene_Name), file = "ensembl_ID_background.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)
write.table(gene_count_h$SE.Gene_Name, file = "high_ensembl_ID.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)

#Export for paper
#high_hgnc <- as.data.frame(high_hgnc)
gene_count_h <- as.data.frame(gene_count_h)
gene_count_h <- gene_count_h[order(gene_count_h$n),]
gt_high <- gt(gene_count_h)
gt_high <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 2.** Ensembl gene IDs containing >5 genetic burden variants")
  ) %>%
  cols_label(SE.Gene_Name = md("**Ensembl gene ID**"),
             n = md("**Number**")) %>%
  cols_width(
    SE.Gene_Name ~ px(100),
    n ~ px(50)
  ) %>% tab_options(
    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(200)
  )
gt_high %>% gtsave("Figures_tables/supp_table2.pdf")

#Get genes with AF >50%
#Get AF in each gene
gene_count2 <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  summarize(AF = mean(AF))
mean(gene_count2$AF)
range(gene_count2$AF)

gene_count2_h <- gene_count2 %>%
  filter(AF > 0.05)

gene_count2_h <- as.data.frame(gene_count2_h)
gene_count2_h <- gene_count2_h[order(gene_count2_h$AF),]
length(unique(gene_count2_h$SE.Gene_Name))
#high_AF_hgnc <- hgnc_count2 %>%
#  filter(AF >0.05)
#length(high_AF_hgnc$AF)
write.table(gene_count2_h$SE.Gene_Name, file = "high_AF_Ensembl.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)

intersect(gene_count_h$SE.Gene_Name, gene_count2_h$SE.Gene_Name)

#Export for DAVID
#high_hgnc$hgnc_symbol
#high_hgnc <- merge(high_hgnc, SE.gene.ID, by = "hgnc_symbol")
#high_AF_hgnc <- merge(high_AF_hgnc, SE.gene.ID, by = "hgnc_symbol")


length(intersect(gene_count2_h$SE.Gene_Name, gene_count_h$SE.Gene_Name))

#Export for paper
gt_high_AF <- gt(gene_count2_h)
gt_high_AF <- 
  gt_high_AF %>%
  tab_header(
    title = md("**Supplementary table 3.** Ensembl gene IDs containing a genetic burden variants at an allele frequency >5%")
  ) %>%
  cols_label(SE.Gene_Name = md("**Ensembl gene ID**"),
             AF = md("**Allele frequency**")) %>%
  fmt_number(columns = AF,
             decimals = 2) %>%
  cols_width(
    SE.Gene_Name ~ px(200),
    AF ~ px(100)
  )%>% tab_options(
    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(200)) 
gt_high_AF %>% gtsave("Figures_tables/supp_table3.pdf")

##Frequency of loss of function variants
table(int_all$SE.LOF)
table(int_all$VEP.Consequence)

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
  mutate(across('SE.Annotation', str_replace_all, "_", " ")) %>%
  mutate(across('VEP.Consequence', str_replace_all, "_", " ")) %>%
  mutate(across('SE.Annotation', str_replace, "bidirectional ", "")) %>%
  mutate(across('VEP.Consequence', str_replace, "bidirectional", ""))

table(high_int_tidy$VEP.Consequence)

high_int_tidy <- high_int_tidy %>%
  mutate(consequence_agreement = if_else(SE.Annotation == VEP.Consequence, "yes", "no"))

table(high_int_tidy$consequence_agreement)
table(high_int_tidy$intersect)

high_int_tidy$consequence <- paste(high_int_tidy$SE.Annotation, ":", high_int_tidy$VEP.Consequence)

int_tbl <- high_int_tidy %>%
  group_by(consequence) %>%
  summarize(count=n())

int_tbl[c('SnpEff', 'Ensembl')] <- str_split_fixed(int_tbl$consequence, ' : ', 2)
int_tbl$consequence <- NULL
int_tbl <- int_tbl %>%
  mutate(consequence_agreement = if_else(SnpEff == Ensembl, "yes", "no"))
 
int_tbl <- int_tbl[c('SnpEff', 'Ensembl', 'count', 'consequence_agreement')]
int_tbl <- int_tbl %>%
  arrange(count)

gt_int <- gt(int_tbl)

gt_int <- gt_int %>%
  tab_style(
    style = cell_fill(color = "lightgray"),
    locations = cells_body(
      rows = consequence_agreement == "yes")) %>%
  tab_header(
    title = md("**Table 2. Overlap between SnpEff and Ensembl-VEP consequence.** Predicted consequence reported by SnpEff and Ensembl-VEP for each genetic burden variant")
  ) %>%
  tab_footnote(
    footnote = "Grey shading represents agreement between SnpEff and Ensembl-VEP.") %>%
  cols_label(SnpEff = md("**SnpEff**"),
             Ensembl = md("**Ensembl-VEP**"),
             count = md("**Count**"),
             consequence_agreement = md("consequence_ag")
             ) %>%
  fmt_number(columns = count, sep_mark = ",", use_seps = TRUE, decimals = 0) %>%
  cols_hide(columns = consequence_agreement) %>%  # Hide the specific column
  cols_width(
    SnpEff ~ px(100),
    Ensembl ~ px(100),
    count ~ px(50)
  )  %>% tab_options(
    table.font.size = 14,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(300)) 

gt_int %>% gtsave("Figures_tables/table2.png")
23026/(length(high_int_tidy$consequence_agreement))

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
    title = md("**Table 2.** Overlap between SnpEff and Ensembl-VEP predicted variant consequence")
  )  %>%
  cols_label(consequence = md("**Consequence**"),
             SE_c = md("**SnpEff**"),
             VEP_c = md("**VEP**")) %>%
  fmt_number(columns = 2:3, sep_mark = ",", use_seps = TRUE, decimals = 0) %>%
  cols_width(
    consequence ~ px(150),
    SE_c ~ px(100),
    VEP_c ~ px(100)
  )
gt_SE_VEP %>% gtsave("Figures_tables/table2.png")


###No gene enrichment using enricher
#test <- biomartr::organismBM(organism = "Equus caballus")
#ecaballus_attributes <- 
#  biomartr::organismAttributes("Equus caballus") %>% 
#  filter(dataset == "ecaballus_gene_ensembl")

#attributes_to_retrieve = c("hgnc_id", "entrezgene_id", "go_id", "name_1006", "namespace_1003")

#result_BM <- biomartr::biomart( genes      = high_AF_hgnc$SE.Gene_Name,                  # genes were retrieved using biomartr::getGenome()
#                                mart       = "ENSEMBL_MART_ENSEMBL",                     # marts were selected with biomartr::getMarts()
#                                dataset    = "ecaballus_gene_ensembl",               # datasets were selected with biomartr::getDatasets()
#                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
#                                filters =   "ensembl_gene_id" )# query key
#head(result_BM)  
#all_EC_annotated <- biomartr::biomart(genes = high_int$SE.Gene_Name,
#                                      mart       = "ENSEMBL_MART_ENSEMBL", 
#                                      dataset    = "ecaballus_gene_ensembl",  
#                                      attributes = attributes_to_retrieve,
#                                      filters =  "ensembl_gene_id" )

#EC_t2g <- as.data.frame(result_BM$go_id)
#EC_t2g$go_di <- result_BM$name_1006

#ora_analysis_high <- enricher(gene = result_BM$name_1006, 
#                              universe = all_EC_annotated$name_1006, 
#                              pAdjustMethod = "bonferroni",
#                              minGSSize = 2,
#                              maxGSSize = 500,
#                              TERM2GENE = EC_t2g,
#                              pvalueCutoff = 0.05)

#dotplot(ora_analysis_high)

#ora_analysis_high <- as.data.frame(ora_analysis_high)
#ora_analysis_bp <- pairwise_termsim(ora_analysis_high, method = "JC")
#emapplot(ora_analysis_high, color = "qvalue")

####Get LOF variants for SnpEff and Ensembl VEP ####
##Need to figure it out manually for VEP
table(high_int_tidy$SE.LOF)

high_int_tidy <- high_int_tidy %>%
  mutate(VEP.LOF = str_detect(VEP.Consequence, "frameshift|splice acceptor|
                              splice donor|start lost|stop lost|stop gained"))
cons <- c("Not_LOF", "LOF")
high_int_tidy <- high_int_tidy %>%
  mutate(across("VEP.LOF", str_replace, "TRUE", "yes")) %>%
  mutate(across("VEP.LOF", str_replace, "FALSE", "no"))

high_int_tidy$LOF <- paste(high_int_tidy$SE.LOF, high_int_tidy$VEP.LOF, sep = ":")
table(high_int_tidy$LOF)

#####LOF per individual
ind$chrom_pos <- paste(ind$chrom, ind$pos, sep = ":")
high_int_tidy$chrom_pos <- paste(high_int_tidy$chrom, high_int_tidy$pos, sep = ":")

LOF <- high_int_tidy$chrom_pos[high_int_tidy$LOF == "yes:yes"]
summary(high_int_tidy$AF[high_int_tidy$LOF == "yes:yes"])
high_int <- high_int %>%
  mutate(VEP.LOF = str_detect(VEP.Consequence, "frameshift|splice_acceptor_variant|
                              splice_donor_variant|splice_region|stop_gained"))
cons <- c("Not_LOF", "LOF")
high_int <- high_int %>%
  mutate(across("VEP.LOF", str_replace, "TRUE", "yes")) %>%
  mutate(across("VEP.LOF", str_replace, "FALSE", "no"))

high_int$LOF <- paste(high_int$SE.LOF, high_int$VEP.LOF, sep = ":")
table(high_int$LOF)

#####LOF per individual
ind_details$chrom_pos <- paste(ind_details$CHROM, ind_details$POS, sep = ":")
high_int$chrom_pos <- paste(high_int$chrom, high_int$pos, sep = ":")

LOF <- high_int$chrom_pos[high_int$LOF == "yes:yes"]
mean(high_int$AF[high_int$LOF == "yes:yes"])
range(high_int$AF[high_int$LOF == "yes:yes"])

ind$chrom_pos <- ind_details$chrom_pos

LOF_ind <- ind %>% 
  filter(ind$chrom_pos %in% LOF)
LOF_ind$chrom_pos <- NULL

LOF_ind_b <- LOF_ind %>%
  mutate_all(., str_replace_all, "/", "") %>%
  mutate_all(., str_replace_all, "\\|", "") %>%
  mutate_all(., str_replace, "01", "het") %>%
  mutate_all(., str_replace, "02", "het") %>%
  mutate_all(., str_replace, "12", "het") %>%
  mutate_all(., str_replace, "12", "het") %>%
  mutate_all(., str_replace, "03", "het") %>%
  mutate_all(., str_replace, "03", "het") %>%
  mutate_all(., str_replace, "10", "het") %>%
  mutate_all(., str_replace, "10", "het") %>%
  mutate_all(., str_replace, "13", "het") %>%
  mutate_all(., str_replace, "13", "het") %>%
  mutate_all(., str_replace, "04", "het") %>%
  mutate_all(., str_replace, "04", "het") %>%
  mutate_all(., str_replace, "11", "hom") %>%
  mutate_all(., str_replace, "11", "hom") %>%
  mutate_all(., str_replace, "22", "hom") %>%
  mutate_all(., str_replace, "22", "hom") %>%
  mutate_all(., str_replace, "33", "hom") %>%
  mutate_all(., str_replace, "33", "hom") %>%
  mutate_all(., str_replace, "00", "hom_WT")  

#Get het/hom by individual

het <- colSums(LOF_ind_b == "het")
hom <- colSums(LOF_ind_b == "hom")
summary(hom)

total_variant <- het + hom
summary(total_variant)
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
                     values = high_int_tidy$SE.Gene_Name[high_int_tidy$LOF == "yes:yes"], 
                     mart = ensembl)
VEP.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                      values = high_int_tidy$VEP.Gene[high_int_tidy$LOF == "yes:yes"], 
                     values = high_int$SE.Gene_Name[high_int$LOF == "yes:yes"], 
                     mart = ensembl)
VEP.gene.ID <- getBM( attributes = c("hgnc_symbol", "ensembl_gene_id"),
                      values = high_int$VEP.Gene[high_int$LOF == "yes:yes"], 
                      mart = ensembl)

SE.HGNC <- SE.gene.ID %>%
  subset(nchar(hgnc_symbol) >1)

#Merge with high_int
colnames(SE.gene.ID) <- c("hgnc_symbol", "SE.Gene_Name")
high_int_m <- merge(high_int_tidy[high_int_tidy$LOF == "yes:yes",], 
                    SE.gene.ID, by = "SE.Gene_Name")
high_int_m <- merge(high_int[high_int$LOF == "yes:yes",], SE.gene.ID, by = "SE.Gene_Name")

gene_count <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  count(SE.Gene_Name)

length(unique(gene_count$SE.Gene_Name))
summary(gene_count$n)
table(gene_count$n)
table(gene_count$n > 5)
mean(gene_count$n)
range(gene_count$n)
table(gene_count$n)

hgnc_count <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  count(hgnc_symbol)

#Remove ""
hgnc_count <- hgnc_count[-1,]
length(hgnc_count$hgnc_symbol)

summary(hgnc_count$n)

mean(hgnc_count$n)
range(hgnc_count$n)
table(hgnc_count$n)
table(hgnc_count$n > 5)

#Get AF in each gene
gene_count2 <- high_int_m %>%
  group_by(SE.Gene_Name) %>%
  summarize(AF = mean(AF))


summary(gene_count2$AF)

mean(gene_count2$AF)
range(gene_count2$AF)

hgnc_count2 <- high_int_m %>%
  group_by(hgnc_symbol) %>%
  summarize(AF = mean(AF))


summary(hgnc_count2$AF)

mean(hgnc_count2$AF)
range(hgnc_count2$AF)


#Get genes containing 5 or more variants
high_LOF <- gene_count %>% 
  filter(n >5)
length(unique(high_LOF$SE.Gene_Name))

#Export for DAVID analysis
write.table(unique(gene_count$SE.Gene_Name), file = "LOF_background.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)
write.table(high_LOF$SE.Gene_Name, file = "LOF_high_Ensembl.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)


#Export for paper
high_LOF <- as.data.frame(high_LOF)
gt_high <- gt(high_LOF)
gt_high <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 4.** Ensembl gene IDs containing >5 LOF variants")
  ) %>%
  cols_label(SE.Gene_Name = md("**Ensembl gene ID**"),
             n = md("**Number**")) %>%
  cols_width(

    SE.Gene_Name ~ px(100),
    n ~ px(50)
  )%>% tab_options(
    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(200)) 


gt_high %>% gtsave("Figures_tables/supp_table4.pdf")

#Get genes with AF >50%
high_AF <- gene_count2 %>%
  filter(AF >0.05)
length(unique(high_AF$SE.Gene_Name))
write.table(high_AF$SE.Gene_Name, file = "LOF_high_AF.txt", quote = FALSE, sep = "\t", row.names = F, col.names = F)

#Export for DAVID

#high_hgnc$hgnc_symbol
#high_hgnc <- merge(high_hgnc, SE.gene.ID, by = "hgnc_symbol")
#high_AF_hgnc <- merge(high_AF_hgnc, SE.gene.ID, by = "hgnc_symbol")

high_hgnc$hgnc_symbol
high_hgnc <- merge(high_hgnc, SE.gene.ID, by = "hgnc_symbol")
high_AF_hgnc <- merge(high_AF_hgnc, SE.gene.ID, by = "hgnc_symbol")


#Export for paper
high_AF <- as.data.frame(high_AF)
gt_high_AF <- gt(high_AF)
gt_high_AF <- 
  gt_high %>%
  tab_header(
    title = md("**Supplementary table 5.** Ensembl gene IDs containing LOF variants at an allele frequency >5%")
  ) %>%
  cols_label(SE.Gene_Name = md("**Gene ID**"),
             n = md("**Number**")) %>%
  cols_width(

    SE.Gene_Name ~ px(100),
    n ~ px(50)
  )%>% tab_options(
    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL,
    table.width = px(200)) 

gt_high_AF %>% gtsave("Figures_tables/supp_table5.pdf")


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

hand2 <- hand %>%
  filter(hand$chrom_pos %in% intersect(hand$chrom_pos, high_int_tidy$chrom_pos))




high_int_tidy$chrom_pos <- paste(high_int_tidy$chrom, high_int_tidy$pos, sep = ":")
high_int_hand <- high_int_tidy %>%
  mutate(hand = if_else(high_int_tidy$chrom_pos == hand$chrom_pos, "yes", "no"))

high_int_hand <- high_int_tidy$chrom_pos[hand$chrom_pos]
high_int_hand <- as.factor(high_int_hand)


####GB and LOF variants by breed####
####Figure out differences in the number of variants per breed
###Individual GB
GB <- ind
GB$chrom_pos <- NULL
ID_GB <- colnames(GB)
het_GB <- data.frame(colSums(GB == "het"))
colnames(het_GB) <- "het_GB"
hom_GB <- as.data.frame(colSums(GB == "hom"))
colnames(hom_GB) <- "hom_GB"
het_LOF <- as.data.frame(colSums(LOF_ind_b == "het"))
colnames(het_LOF) <- "het_LOF"
hom_LOF <- as.data.frame(colSums(LOF_ind_b == "hom"))
colnames(hom_LOF) <- "hom_LOF"

GB.df <- cbind(ID_GB, het_GB$het_GB, hom_GB$hom_GB, het_LOF$het_LOF, hom_LOF$hom_LOF)
GB.df <- as.data.frame(GB.df)
colnames(GB.df) <- c("ID", "het_GB", "hom_GB", "het_LOF", "hom_LOF")
GB.df$het_GB <- as.numeric(GB.df$het_GB)
GB.df$hom_GB <- as.numeric(GB.df$hom_GB)
GB.df$het_LOF <- as.numeric(GB.df$het_LOF)
GB.df$hom_LOF <- as.numeric(GB.df$hom_LOF)

GB.df$GB <- GB.df$het_GB + GB.df$hom_GB
GB.df$LOF <- GB.df$het_LOF + GB.df$hom_LOF


#Add in breed information
breed.df <- bcftools_br[c("ID", "Breed", "DOC")]

GB_b <- merge(GB.df, breed.df, "ID")

#Check for best model
fit1 <- (lm(GB ~ Breed, data=GB_b))
fit2 <- (lm(GB ~ Breed + DOC, data=GB_b))
anova(fit1,fit2) #fit2

#Emmeans - GB
gb_m <- (lm(GB ~ Breed + DOC,data=GB_b))
hom_gb_m <- (lm(hom_GB ~ Breed + DOC,data=GB_b))

#Get EMMEANs
gb_emm <- emmeans(gb_m, specs = "Breed", weights = "proportional", 
                  type = "response")

gb_emm <- as.data.frame(gb_emm)
gb_emm <- gb_emm %>%
  arrange(emmean)

gb_hom_emm <- emmeans(hom_gb_m, specs = "Breed", weights = "proportional", 
                  type = "response")
gb_hom_emm <- as.data.frame(gb_hom_emm)
gb_hom_emm <- gb_hom_emm %>%
  arrange(emmean)

x1 <- plot(gb_emm) + geom_boxplot(col = color_mapping[gb_emm$Breed]) + 
  theme_bw() + xlab("EMMEAN of genetic burden") + 


gb_hom_emm <- emmeans(hom_gb_m, specs = "Breed", weights = "proportional", 
                  type = "response")

x1 <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of genetic burden") + 

  ylab("Breed") +
  scale_x_continuous(breaks = seq(550,1000,50))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

x2 <- plot(gb_hom_emm) + geom_boxplot(col = color_mapping[gb_hom_emm$Breed]) + theme_bw() + 
  xlab("EMMEAN of homozygous genetic burden") + 
  ylab("Breed") +
  scale_x_continuous(breaks = seq(150,400,25))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

gb_t <- gb_emm %>%
  arrange(emmean)

gb_hom_t <- gb_hom_emm %>%
  arrange(emmean)

gb_full <- rbind(gb_t, gb_hom_t)
gb_full <- gt(gb_full)


gb_full <- 
  gb_full %>%
  tab_header(
    title = md("**Table 3. Estimated marginal means (EMMEANs) of the genetic burden by breed.** 
               EMMEAN accounting for depth of coverage, standard error (SE), and 95% confidence intervals (CI) for the genetic 
               burden for each breed with 17 or greater individuals for all variants and those 
               only present in homozygous states.")
  ) %>%
  opt_align_table_header(align = "left") %>%
  cols_label(Breed = md("**Breed**"),
             emmean = md("**EMMEAN**"),
             SE = md("**SE**"),
             lower.CL = md("**Lower CI**"),
             upper.CL = md("**Upper CI**")
             ) %>%
  tab_row_group(
    label = md("**All genetic burden variants**"),
    rows = 1:13) %>% 
  tab_options(
      table.font.size = 8) %>%
  tab_row_group(
    label = md("**Homozygous variants**"),
    rows = 14:26
  ) %>% 
  tab_style(
    style = list(
      cell_fill(color = "gray"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == max(emmean[14:26]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "gray"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == max(emmean[1:13]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "#8aedea"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == min(emmean[1:13]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "cyan"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == min(emmean[14:26]))) %>%
    tab_style(
      style = cell_borders(
        sides = "all", color = "#000000", style = "solid", weight = px(1.5)),
      locations = cells_body(
        columns = everything(),
        rows = everything()
      )
    ) %>%
  fmt_number(columns = c(emmean, SE, lower.CL, upper.CL), sep_mark = ",", use_seps = TRUE, decimals = 0) %>%
  tab_footnote(
    footnote = "Maximum EMMEAN (grey) and minimum EMMEAN (cyan)"
  ) %>%
  cols_hide(
    column = df
  ) %>%
  cols_align(
    align = "left",
    column = Breed) %>%
  tab_options(
    table.font.size = 8
  ) %>%
  cols_width(
    Breed ~ px(90),
    emmean ~px(40),
    SE ~ px(25),
    lower.CL ~ px(40),
    upper.CL ~ px(40)) %>%
  tab_options(data_row.padding = px(5),
              table.font.color = "black",
              table.font.color.light = NULL,
              table.width = px(250))
  
gb_full %>% gtsave("Figures_tables/table3.png")

######LOF table
#Check for best model
fit1 <- (lm(LOF ~ Breed, data=GB_b))
fit2 <- (lm(LOF ~ Breed + DOC, data=GB_b))
anova(fit1,fit2) #fit2

#Emmeans - GB
LOF_m <- (lm(LOF ~ Breed + DOC,data=GB_b))
hom_LOF_m <- (lm(hom_LOF ~ Breed + DOC,data=GB_b))

#Get EMMEANs
lof_emm <- emmeans(LOF_m, specs = "Breed", weights = "proportional", 
                  type = "response")
lof_emm <- as.data.frame(lof_emm)
lof_emm <- lof_emm %>%
  arrange(emmean)

lof_hom_emm <- emmeans(hom_LOF_m, specs = "Breed", weights = "proportional", 
                      type = "response")
lof_hom_emm <- as.data.frame(lof_hom_emm)
lof_hom_emm <- lof_hom_emm %>%
  arrange(emmean)

x3 <- plot(lof_emm) + geom_boxplot(col = color_mapping[lof_emm$Breed]) + 
  theme_bw() + xlab("EMMEAN of LOF variants") + 
  ylab("Breed") +
  scale_x_continuous(breaks = seq(100,250,25))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
x4 <- plot(lof_hom_emm) + geom_boxplot(col = color_mapping[lof_hom_emm$Breed]) + theme_bw() + 
  xlab("EMMEAN of homozygous LOF variants") + 

x2 <- plot(gb_hom_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of homozygous genetic burden") + 
  ylab("Breed") +
  scale_x_continuous(breaks = seq(150,400,50))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
x3 <- plot(lof_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of LOF variants") + 
  ylab("Breed") +
  scale_x_continuous(breaks = seq(250,600,50))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
x4 <- plot(lof_hom_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of homozygous LOF variants") + 

  ylab("Breed") +
  scale_x_continuous(breaks = seq(50,225,25))+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

c <- ggarrange(x1, x2, x3, x4, ncol = 1, nrow = 4, labels = "auto")#, common.legend = TRUE, legend = "bottom")
save_plot("Figures_tables/GB_EMMEANs.pdf", c,
          base_height = 16, base_width = 8)


######LOF table
#Check for best model
fit1 <- (lm(LOF ~ Breed, data=GB_b))
fit2 <- (lm(LOF ~ Breed + DOC, data=GB_b))
anova(fit1,fit2) #fit2

#Emmeans - GB
LOF_m <- (lm(LOF ~ Breed + DOC,data=GB_b))
hom_LOF_m <- (lm(hom_LOF ~ Breed + DOC,data=GB_b))

#Get EMMEANs
lof_emm <- emmeans(LOF_m, specs = "Breed", weights = "proportional", 
                  type = "response")

lof_hom_emm <- emmeans(hom_LOF_m, specs = "Breed", weights = "proportional", 
                      type = "response")
lof_t <- gt(summary(lof_emm))
lof_hom_t <- gt(summary(lof_hom_emm))

lof_full <- rbind(lof_t$`_data`, lof_hom_t$`_data`)
lof_full <- gt(lof_full)

lof_full <- 
  lof_full %>%
  tab_header(
    title = md("**Table 4.** Estimated marginal means (EMMEANs) of the LOF genetic burden by breed"),
    subtitle = md("EMMEAN, standard error (SE), and 95% confidence intervals (CI) for the LOF genetic burden for each breed with 17 or greater individuals for all variants and those only present in homozygous states.")
  ) %>%
  opt_align_table_header(align = "left") %>%
  cols_label(Breed = md("**Breed**"),
             emmean = md("**EMMEAN**"),
             SE = md("**SE**"),
             lower.CL = md("**Lower CI**"),
             upper.CL = md("**Upper CI**")
  ) %>%
  tab_row_group(
    label = md("**All LOF variants**"),
    rows = 1:13) %>%
  tab_row_group(
    label = md("**Homozygous variants**"),
    rows = 13:26
  ) %>% 
  tab_style(
    style = list(
      cell_fill(color = "gray"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == max(emmean[14:26]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "gray"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == max(emmean[1:13]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "cyan"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == min(emmean[1:13]))) %>%
  tab_style(
    style = list(
      cell_fill(color = "cyan"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = emmean,
      rows = emmean == min(emmean[14:26]))) %>%

  tab_style(
    style = cell_borders(
      sides = "all", color = "#000000", style = "solid", weight = px(1.5)),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) %>%
  fmt_number(columns = c(emmean, SE, lower.CL, upper.CL), sep_mark = ",", use_seps = TRUE, decimals = 0) %>%
  tab_footnote(
    footnote = "Maximum EMMEAN (grey) and minimum EMMEAN (cyan)"
  ) %>%
  cols_hide(
    column = df
  ) %>%
  cols_align(
    align = "left",
    column = Breed) %>%
  tab_options(

    table.font.size = 8,
    table.font.color = "black",
    table.font.color.light = NULL
  ) %>%
  cols_width(
    Breed ~ px(90),
    emmean ~px(40),
    SE ~ px(25),
    lower.CL ~ px(40),
    upper.CL ~ px(40)) %>%
  tab_options(data_row.padding = px(5),
              table.width = px(250))

lof_full %>% gtsave("Figures_tables/table4.png")



######Look at effective population size######
# Add in Ne from Jessica's paper
GB_b$ne_JP <- "NA"
GB_b$ne_JP[GB_b$Breed == "Arabian"] <- 346
GB_b$ne_JP[GB_b$Breed == "Belgian"] <- 431
GB_b$ne_JP[GB_b$Breed == "Clydesdale"] <- 194

GB_b$ne_JP[GB_b$Breed == "Franchese Montagne"] <- 316

GB_b$ne_JP[GB_b$Breed == "Franches Montagne"] <- 316

GB_b$ne_JP[GB_b$Breed == "Shetland"] <- 365
GB_b$ne_JP[GB_b$Breed == "Icelandic"] <- 555
GB_b$ne_JP[GB_b$Breed == "Morgan"] <- 448
GB_b$ne_JP[GB_b$Breed == "Quarter Horse"] <- 426
GB_b$ne_JP[GB_b$Breed == "Standardbred"] <- 290
GB_b$ne_JP[GB_b$Breed == "Thoroughbred"] <- 190
GB_b$ne_JP <- as.numeric(GB_b$ne_JP)


JP_GB <- cor.test(x=GB_b$ne_JP,y=GB_b$GB, method = "pearson", use='complete.obs')
JP_hom_GB <- cor.test(x=GB_b$ne_JP,y=GB_b$hom_GB, method = "pearson", use='complete.obs')
JP_LOF <- cor.test(x=GB_b$ne_JP,y=GB_b$LOF, method = "pearson", use='complete.obs')
JP_hom_LOF <- cor.test(x=GB_b$ne_JP,y=GB_b$hom_LOF, method = "pearson", use='complete.obs')


cor.test(x=GB_b$ne_JP,y=GB_b$GB, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_JP,y=GB_b$hom_GB, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_JP,y=GB_b$LOF, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_JP,y=GB_b$hom_LOF, method = "pearson", use='complete.obs')


#Add in effective population size (from Sam's paper)
GB_b$ne_SB <- "NA"
GB_b$ne_SB[GB_b$Breed == "Arabian"] <- 3561
GB_b$ne_SB[GB_b$Breed == "Belgian"] <- 3570
GB_b$ne_SB[GB_b$Breed == "Icelandic"] <- 2736

GB_b$ne_SB[GB_b$Breed == "Franchese Montagne"] <- 3305

GB_b$ne_SB[GB_b$Breed == "Franches Montagne"] <- 3305

GB_b$ne_SB[GB_b$Breed == "Morgan"] <- 4481
GB_b$ne_SB[GB_b$Breed == "Icelandic"] <- 2736
GB_b$ne_SB[GB_b$Breed == "Quarter Horse"] <- 6516
GB_b$ne_SB[GB_b$Breed == "Standardbred"] <- 2528
GB_b$ne_SB[GB_b$Breed == "Thoroughbred"] <- 1784
GB_b$ne_SB[GB_b$Breed == "Welsh Pony"] <- 5625

GB_b$ne_SB <- as.numeric(GB_b$ne_SB)

SB_GB <- cor.test(x=GB_b$ne_SB,y=GB_b$GB, method = "pearson", use='complete.obs')
SB_hom_GB <- cor.test(x=GB_b$ne_SB,y=GB_b$hom_GB, method = "pearson", use='complete.obs')
SB_LOF <- cor.test(x=GB_b$ne_SB,y=GB_b$LOF, method = "pearson", use='complete.obs')
SB_hom_LOF <- cor.test(x=GB_b$ne_SB,y=GB_b$hom_LOF, method = "pearson", use='complete.obs')

#cor.test(x=GB_b$ne_SB,y=GB_b$ne_JP, method = "pearson", use='complete.obs')

corr_JP_SB <- data.frame(
  Ne = c("JP", "JP", "JP", "JP", "SB", "SB", "SB", "SB"),
  Estimate = c("GB", "hom_GB", "LOF", "hom_LOF","GB", "hom_GB", "LOF", "hom_LOF"),
  corr = c(JP_GB$estimate, JP_hom_GB$estimate, JP_LOF$estimate, JP_hom_LOF$estimate, SB_GB$estimate, SB_hom_GB$estimate, SB_LOF$estimate, SB_hom_LOF$estimate),
  low_conf = c(JP_GB$conf.int[1], JP_hom_GB$conf.int[1], JP_LOF$conf.int[1], JP_hom_LOF$conf.int[1],SB_GB$conf.int[1], SB_hom_GB$conf.int[1], SB_LOF$conf.int[1], SB_hom_LOF$conf.int[1]),
  high_conf = c(JP_GB$conf.int[2], JP_hom_GB$conf.int[2], JP_LOF$conf.int[2], JP_hom_LOF$conf.int[2], SB_GB$conf.int[2], SB_hom_GB$conf.int[2], SB_LOF$conf.int[2], SB_hom_LOF$conf.int[2]),
  p = c(JP_GB$p.value, JP_hom_GB$p.value, JP_LOF$p.value, JP_hom_LOF$p.value, SB_GB$p.value, SB_hom_GB$p.value, SB_LOF$p.value, SB_hom_LOF$p.value)
)


x1 <- ggplot(corr_JP_SB, aes(x = Estimate, y = corr, color = Ne)) + 
  geom_point(position = position_dodge(width = 0.8), size = 3) + 
  geom_errorbar(aes(ymin = low_conf, ymax = high_conf),
                position = position_dodge(width = 0.8), size = 0.3) +
  theme_bw() + xlab("Variant type") + 
  ylab("Pearson's correlation coeeficient") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

save_plot("Figures_tables/NE_results.png", x1,
          base_height = 8, base_width = 8)

########################Known variants

cor.test(x=GB_b$ne_SB,y=GB_b$GB, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_SB,y=GB_b$hom_GB, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_SB,y=GB_b$LOF, method = "pearson", use='complete.obs')
cor.test(x=GB_b$ne_SB,y=GB_b$hom_LOF, method = "pearson", use='complete.obs')



#########OMIA variants############# Didn't do this because needed to check the exact position of the variant
#Variants <= 20 bp in length
snpMart = useEnsembl(biomart = "snps",
                     dataset = "ecaballus_snp")

##Known variants

##Have to tidy the files outside of excel so it doesn't mess things up
known_gt <- read.table("known_variants_tidy.txt", header = T)
known_pos <- read.table("known_variants.txt", header = T)

#Check variants are the same
known_gt$chrom_pos <- paste(known_gt$chrom, known_gt$pos, sep = ":")
known_pos$chrom_pos <- paste(known_pos$Chromosome, known_pos$Start, sep = ":")
known_gt$M989 <- NULL
known_gt$M6468 <- NULL

known_pos$overlap <- known_gt$chrom_pos == known_pos$chrom_pos
known_gt$chrom_pos[known_pos$overlap == "FALSE"]

#Need to make this into a table with total no. OMIA variants,
#total no in this pop, no. heterozygotes, no. homozygotes, 
#mean (range) AF
known_pos$type <- paste(known_pos$Causative., known_pos$Deleterious., sep = ":")

table(known_pos$type) #causative:deleterious
colnames(known_pos)[1] <- "phenotype"

known_all <- merge(known_pos, known_gt, "phenotype")
summary(known_all$AF, na.rm = T)
table(known_all$type[known_all$chrom != "NA"])

known_tbl1 <- known_all %>%
  select(phenotype, AF, Chromosome, type)

known_tbl2 <- known_tbl1 %>%
  group_by(type) %>%
  summarize(count=n(), AF=median(AF, na.rm = T)*100)

OMIA <- as.data.frame(table(known_all$type[known_all$chrom != "NA"]))
colnames(OMIA) <- c("type", "OMIA_count")
IQR1 <- known_tbl1 %>%
  group_by(type) %>%
  summarize(IQR1=summary(AF, na.rm = T)[2]*100)
IQR1<- as.data.frame(IQR1)
colnames(IQR1) <- c("type", "IQR1")

IQR3 <- known_tbl1 %>%
  group_by(type) %>%
  summarize(IQR3=summary(AF, na.rm = T)[5]*100)
IQR3<- as.data.frame(IQR3)
colnames(IQR3) <- c("type", "IQR3")

known_tbl1 <- merge(known_tbl2, OMIA, "type")
known_tbl1 <- merge(known_tbl1, IQR1, "type") 
known_tbl1 <- merge(known_tbl1, IQR3, "type") 

known_tbl1 <- known_tbl1 %>%
  select(type, count, OMIA_count, AF, IQR1, IQR3)

x <- c("yes:yes", "CC:CC", "no:yes", "no:no")
known_tbl3 <- known_tbl1 %>%
  mutate(type =  factor(type, levels = x)) %>%
  arrange(type)

known_tbl4 <- known_tbl3 %>%
  mutate(across('type', str_replace, "no:no", "Non-disease and non-coat color trait associated")) %>%
  mutate(across('type', str_replace, "CC:CC", "Coat color causing and associated")) %>%
  mutate(across('type', str_replace, "no:yes", "Disease associated")) %>%
  mutate(across('type', str_replace, "yes:yes", "Disease and non-coat color trait causing"))
known_tbl4$IQR1 <- as.numeric(known_tbl4$IQR1)
known_tbl4$IQR3 <- as.numeric(known_tbl4$IQR3)
known_tb_gt <- gt(known_tbl4)

known_tb_gt <- 
  known_tb_gt %>%
  fmt_number(columns = c(AF, IQR1, IQR3), sep_mark = ",", use_seps = TRUE, decimals = 2) %>%
  cols_merge_range(IQR1, IQR3, sep = " - ") %>% 
  tab_header(
    title = md("**Table 5. Classification of OMIA variants by type.** Number of: causative variants for disease and non-coat color traits, associated and causative variants for coat color, associated variants for disease, and associated variants for non-disease and non-coat color traits present in this population, and the median and range of the allele frequency.")
  ) %>%
  opt_align_table_header(align = "left") %>%
  cols_label(type = md("**OMIA variant type**"),
             count = md("**Number of OMIA variants**"),
             OMIA_count = md("**Number of OMIA variants in this population**"),
             AF = md("**Median allele frequency (%)**"),
             IQR1 = md("**AF interquartile range (%)**")
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "all", color = "#000000", style = "solid", weight = px(1.5)),
    locations = cells_body(
      columns = everything(),
      rows = everything()
    )
  ) %>%
  cols_align(
    align = "left",
    column = type) %>%
    cols_align(
      align = "center",
      column = count) %>%
    cols_align(
      align = "center",
      column = OMIA_count) %>%
    cols_align(
      align = "center",
      column = AF) %>%
    cols_align(
      align = "center",
      column = IQR1) %>%
  tab_options(
    table.font.size = 12
  ) %>%
  cols_width(
    type ~ px(110),
    count ~px(100),
    OMIA_count ~ px(100),
    AF ~ px(100),
    IQR1 ~ px(100)) %>%
  tab_options(data_row.padding = px(5))%>% 
  tab_options(table.font.size = 10,
              table.font.color = "black",
              table.font.color.light = NULL) 

known_tb_gt %>% gtsave("Figures_tables/table5.png", expand = 0)


#Get figure 3 - known variants present
known_s <- read.table("known_variants_summary_by_breed.txt", header = T)
colnames(known_s) <- c("phenotype", "Breed", "genotype", "genotype_count", "AC", "AF")
known_sp <- merge(known_s, known_pos, "phenotype")

#Causative:deleterious
#Plot
known_yy <- known_sp %>%
  filter(type == "yes:yes") %>%
  filter(genotype_count != 0)

yy <- ggplot(known_yy, aes(x = phenotype, y = genotype_count, 
                           fill = Breed, group = genotype)) + 
  geom_bar_pattern(stat = "identity", aes(pattern = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) + 

  scale_fill_manual(values = c("grey", "blue", "purple", "orchid1", "khaki", "dodgerblue",
                               "green4", "black", "maroon", "#FDBF6F", "yellow")) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,50)) + 
  xlab("Disease and non-coat color trait causing variants") + 

  scale_fill_manual(values = as.vector(alphabet2(13)),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,25)) + 
  xlab("Disease causing variants") + 

  theme(axis.text = element_text(size=12,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=14,face="bold"))


known_CC <- known_sp %>%
  filter(type == "CC:CC") %>%
  filter(genotype_count != 0)

CC <- ggplot(known_CC, aes(x = phenotype, y = genotype_count,
                           fill = Breed, group = genotype)) +
  geom_bar_pattern(stat = "identity", aes(pattern = factor(genotype)),
                 position = position_dodge(width = 1)) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) + 

  scale_fill_manual(values = c("grey", "blue", "cyan", "red", "purple", "orchid1", 
                               "khaki", "dodgerblue", "green4", "black", "maroon",
                               "#FDBF6F", "yellow")) +
  scale_alpha_manual(values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,100)) + 
  xlab("Coat color associated and causative variants") + 

  scale_fill_manual(values = as.vector(alphabet2(13)),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_alpha_manual(values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,125)) + 
  xlab("Non-disease causing variants") + 

  theme(axis.text = element_text(size=12,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=14,face="bold"))

known_ny <- known_sp %>%
  filter(type == "no:yes") %>%
  filter(genotype_count != 0)

ny <- ggplot(known_ny, aes(x = phenotype, y = genotype_count, 
                           fill = Breed, group = genotype)) + 
  geom_bar_pattern(stat = "identity", aes(pattern = factor(genotype)),
                   position = position_dodge(width = 1)) +
  scale_pattern_manual(values = c("none", "stripe")) + 

  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) +  
  scale_fill_manual(values = c("grey", "blue", "cyan", "red", "purple", "orchid1", 
                               "khaki", "dodgerblue", "green4", "black", "maroon",
                               "#FDBF6F", "yellow")) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,100)) + 

  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) + 
  scale_fill_manual(values = as.vector(alphabet2(13)), 
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,125)) + 

  xlab("Disease associated variants") + 
  theme(axis.text = element_text(size=12,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=14,face="bold"))

known_nn <- known_sp %>%
  filter(type == "no:no") %>%
  filter(genotype_count != 0)

nn <- ggplot(known_nn, aes(x = phenotype, y = genotype_count, 
                           fill = Breed, group = genotype)) + 
  geom_bar_pattern(stat = "identity", aes(pattern = factor(genotype)),
                   position = position_dodge(width = 1)) +
  scale_pattern_manual(values = c("none", "stripe")) + 

  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) +  
  scale_fill_manual(values = c("grey", "blue", "cyan", "red", "purple", "orchid1", 
                               "khaki", "dodgerblue", "green4", "black", "maroon",
                               "#FDBF6F", "yellow"), 
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,80)) + 
  xlab("Non-disease and non-coat color trait associated variants") + 

  scale_pattern_density_manual(values = c(het = 0, hom = 0.05)) + 
  scale_fill_manual(values = as.vector(alphabet2(13)),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,80)) + 
  xlab("Non-disease associated variants") + 

  theme(axis.text = element_text(size=12,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size=14,face="bold"))

#Save as combined plot

OMIA_fig <- ggarrange(yy, CC, ny, nn, ncol=1, nrow=4, common.legend = TRUE, 
                  legend = "bottom", labels = "auto")
save_plot("Figures_tables/OMIA_variants.png", 
          OMIA_fig, base_height = 16,base_width = 12)

#Get the legend from another figure
nn_p <- ggarrange(nn, ncol=1, nrow=1, legend = "bottom")
save_plot("Figures_tables/OMIA_legend.png", nn_p, base_height = 6, base_width = 12)



#Get table of allele frequencies by breed
known_sp_tidy <- known_sp[!duplicated(known_sp[c(1,2)]),]
known_sp_tidy <- known_sp_tidy %>% arrange(Breed)

#Get the OMIA ID and full phenotype name
details <- read.table("known_variants_details.txt", header = TRUE)
known_sp_tidy2 <- merge(known_sp_tidy, details, col = "phenotype")

known_sp_tidy3 <- known_sp_tidy2[,c("Phenotype_full", "phenotype", "Breed", "AF", 
                                   "Gene", "Deleterious.", "Causative.", 
                                   "Ref", "Alt", "chrom_pos")]
known_sp_tidy3$Phenotype_full <- gsub("_", " ", known_sp_tidy3$Phenotype_full, fixed=TRUE)

gt_known_sp <- gt(known_sp_tidy3)
gt_known_sp <- 
  gt_known_sp %>%
  tab_header(
    title = md("**Supplementary table 6.** Breed specific allele frequencies for the known variants."),
    subtitle = "Variants not present in this population are listed as NA for the breed and allele frequency"
  ) %>%
  cols_label(
             Phenotype_full = md("**Phenotype***"),
             Breed = md("**Breed**"),
             phenotype = md("**Phenotype abbreviation***"),
             AF = md("**Allele frequency**"),
             Gene = md("**Gene name**"),
             Deleterious. = md("**Disease**"),
             Causative. = md("**Causative**"),
             Ref = md("**Reference allele**"),
             Alt = md("**Alternate allele**"),
             chrom_pos = md("**chromosome:position**")) %>%
  fmt_number(columns = AF, sep_mark = ",", use_seps = TRUE, decimals = 3) %>%
  cols_align(align = "center",
             columns = everything()) %>%
  cols_width(
    Phenotype_full ~ px(100),
    phenotype ~ px(100),
    Breed ~ px(75),
    AF ~ px(50),
    Gene ~ px(50),
    Deleterious. ~ px(75),
    Causative. ~ px(75),
    Ref ~ px(50),
    Alt ~ px(50),
    chrom_pos ~ px(100)) %>%
  tab_source_note(
    source_note = "NA: Variants not present in this population") %>%
  tab_source_note(
    source_note = "*Phenotype and abbreviation pulled from OMIA"
  ) %>%
  tab_source_note(
    source_note = "CC = coat color related"
  ) %>% 
  tab_options(
    table.font.size = 8) 
gt_known_sp %>% gtsave("Figures_tables/supp_table6.pdf", expand = 0)


###Get table of horse IDs for IMM, GBED, gait, and CA
##IMM
rownames(known_all) <- known_all$phenotype
known_IMM <- known_all %>%
  filter(rownames(known_all) == "MYHM") %>%
  select(where(~ grepl("0/1",.)))

known_IMM_b <- t(known_IMM) 
known_IMM_b <- as.data.frame(known_IMM_b)
known_IMM_b$ID <- rownames(known_IMM_b)
known_IMM_b <- merge(known_IMM_b, breeds, "ID")

##CLF
rownames(known_all) <- known_all$phenotype
known_CLF <- known_all[11:15,]

selected_columns <- known_CLF[, apply(known_CLF, 2, function(x) any(grepl("0/1", x) | grepl("1/1", x)))]

selected_columns2 <- selected_columns %>%
 mutate_all(., str_replace_all, "/", "") %>%
  mutate_all(., str_replace, "01", "1") %>%
  mutate_all(., str_replace, "11", "1")

to_numeric <- function(x) {
  as.numeric(as.character(x))
}
# Apply the function to all columns in the data frame
data_numeric <- apply(selected_columns2, 2, to_numeric)

# Convert the resulting matrix back to a data frame
data_numeric <- as.data.frame(data_numeric)

sel_cols <- colSums(data_numeric[1:369])
which(sel_cols ==5) 

known_CLF %>%
  filter(rownames(known_all) == "CLF") %>%
  select(where(~ grepl("0/1",.)))

known_IMM_b <- t(known_IMM) 
known_IMM_b <- as.data.frame(known_IMM_b)
known_IMM_b$ID <- rownames(known_IMM_b)
known_IMM_b <- merge(known_IMM_b, breeds, "ID")

##GBED
known_GBED <- known_all %>%
  filter(rownames(known_all) == "GBED") %>%
  select(where(~ grepl("1/1",.)))

known_GBED_b <- t(known_GBED) 
known_GBED_b <- as.data.frame(known_GBED_b)
known_GBED_b$ID <- rownames(known_GBED_b)
known_GBED_b <- merge(known_GBED_b, breeds, "ID")

##CA
known_CA <- known_all %>%
  filter(rownames(known_all) == "CA") %>%
  select(where(~ grepl("0/1",.)))

known_CA_b <- t(known_CA) 
known_CA_b <- as.data.frame(known_CA_b)
known_CA_b$ID <- rownames(known_CA_b)
known_CA_b <- merge(known_CA_b, breeds, "ID")

#Gait
known_gait <- known_all %>%
  filter(rownames(known_all) == "gait") %>%
  select(where(~ grepl("0/1",.)))

known_gait_b <- t(known_gait) 
known_gait_b <- as.data.frame(known_gait_b)
known_gait_b$ID <- rownames(known_gait_b)
known_gait_b <- merge(known_gait_b, breeds, "ID")

#hydro
known_hydro <- known_all %>%
  filter(rownames(known_all) == "hydro") %>%
  select(where(~ grepl("0/1",.)))

known_hydro_b <- t(known_hydro) 
known_hydro_b <- as.data.frame(known_hydro_b)
known_hydro_b$ID <- rownames(known_hydro_b)
known_hydro_b <- merge(known_hydro_b, breeds, "ID")

#Coat sun coat_sun
known_coat_sun <- known_all %>%
  filter(rownames(known_all) == "coat_sun") %>%
  select(where(~ grepl("1/1",.)))

known_coat_sun_b <- t(known_coat_sun) 
known_coat_sun_b <- as.data.frame(known_coat_sun_b)
known_coat_sun_b$ID <- rownames(known_coat_sun_b)
known_coat_sun_b <- merge(known_coat_sun_b, breeds, "ID")
length(known_coat_sun_b$coat_sun)

#PSSM1 - not interesting
known_PSSM1 <- known_all %>%
  filter(rownames(known_all) == "PSSM1") %>%
  select(where(~ grepl("1/1",.)))

known_PSSM1_b <- t(known_PSSM1) 
known_PSSM1_b <- as.data.frame(known_PSSM1_b)
known_PSSM1_b$ID <- rownames(known_PSSM1_b)
known_PSSM1_b <- merge(known_PSSM1_b, breeds, "ID")

##FFS
known_FFS <- known_all %>%
  filter(rownames(known_all) == "FFS") %>%
  select(where(~ grepl("0/1",.)))

known_FFS_b <- t(known_FFS) 
known_FFS_b <- as.data.frame(known_FFS_b)
known_FFS_b$ID <- rownames(known_FFS_b)
known_FFS_b <- merge(known_FFS_b, breeds, "ID")

##D2
known_D2 <- known_all %>%
  filter(rownames(known_all) == "dwarf_D2") %>%
  select(where(~ grepl("0/1",.)))

known_D2_b <- t(known_D2) 
known_D2_b <- as.data.frame(known_D2_b)
known_D2_b$ID <- rownames(known_D2_b)
known_D2_b <- merge(known_D2_b, breeds, "ID")

##agouti
known_coat_bl <- known_all %>%
  filter(rownames(known_all) == "coat_bl") %>%
  select(where(~ grepl("1/1",.))) 

known_coat_bl2 <- t(known_coat_bl) 
known_coat_bl2 <- as.data.frame(known_coat_bl2)
known_coat_bl2$ID <- rownames(known_coat_bl2)
known_coat_bl2 <- merge(known_coat_bl2, breeds, "ID")
known_coat_bl2[known_coat_bl2$Breed == "Belgian",]

##chestnut
known_coat_chest12 <- known_all %>%
  filter(rownames(known_all) == "coat_chest1") %>%
  select(where(~ grepl("1/1",.))) 

known_coat_chest12 <- t(known_coat_chest12) 
known_coat_chest12 <- as.data.frame(known_coat_chest12)
known_coat_chest12$ID <- rownames(known_coat_chest12)
known_coat_chest12 <- merge(known_coat_chest12, breeds, "ID")
known_coat_chest12[known_coat_chest12$Breed == "Belgian",]

##brindle
known_coat_bl <- known_all %>%
  filter(rownames(known_all) == "coat_ND1") %>%
  select(where(~ grepl("1/1",.))) 

known_coat_bl2 <- t(known_coat_bl) 
known_coat_bl2 <- as.data.frame(known_coat_bl2)
known_coat_bl2$ID <- rownames(known_coat_bl2)
known_coat_bl2 <- merge(known_coat_bl2, breeds, "ID")
known_coat_bl2

##other coat colors
known_coat_bl <- known_all %>%
  filter(rownames(known_all) == "coat_W32") %>%
  select(where(~ grepl("0/1",.))) 

known_coat_bl2 <- t(known_coat_bl) 
known_coat_bl2 <- as.data.frame(known_coat_bl2)
known_coat_bl2$ID <- rownames(known_coat_bl2)
known_coat_bl2 <- merge(known_coat_bl2, breeds, "ID")
known_coat_bl2


OMIA_fig <- ggarrange(yy, yn, ny, nn, ncol=1, nrow=4, common.legend = TRUE, 
                  legend = "bottom", labels = "auto")
save_plot("Figures_tables/OMIA_variants.pdf", 
          OMIA_fig, base_height = 16,base_width = 12)

#Get the legend from another figure
yy_p <- ggarrange(yy, ncol=1, nrow=1, legend = "bottom")
save_plot("Figures_tables/OMIA_legend.tiff", yy, base_height = 6, base_width = 12)




# Pull out CLF variants for zoom
CLF <- data1 %>%
  filter(grepl(c('CLF'), Phenotype))
bp4 <- ggplot(CLF, aes(x = Phenotype, y = count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,80)) + 
  xlab("Congenital Liver Fibrosis alleles") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"))

mean(CLF$AF)
range(CLF$AF)
