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
library(reshape2)
library(ggrepel)
library(Gviz)
library(devtools)
library(nlcor)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

#Only include autosomes and chr X (not MT and unplaced contigs)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")

####DOC and other descriptive stats
intersect_stats <- read.table("intersect_by_ind_number_of_variants.txt",header=T)
intersect_stats$nvariants <- intersect_stats$nNonRefHom + intersect_stats$nHets
mean(intersect_stats$nvariants)
mean(intersect_stats$nIndels)
mean(intersect_stats$nNonRefHom)
mean(intersect_stats$nHets)
intersect_stats$tstv <- intersect_stats$Ts/intersect_stats$Tv

#Get number of heterozygous variants and number of homozygous variants per kb of sequence
#taken from NCBI (2,474.93 Mb, 2,474,930 kb)
intersect_stats$het_kb <- intersect_stats$nHets /2474930
intersect_stats$hom_kb <- intersect_stats$nNonRefHom /2474930
intersect_stats$nvariants_kb <- intersect_stats$nvariants /2474930
mean(intersect_stats$het_kb)
range(intersect_stats$het_kb)
mean(intersect_stats$hom_kb)
range(intersect_stats$hom_kb)
mean(intersect_stats$nvariants_kb)
range(intersect_stats$nvariants_kb)

#Add in DOC information
DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) = c("Sample", "total_DOC","nuclear_placed_DOC")
intersect_doc <- merge(intersect_stats,DOC, by="Sample")
summary(intersect_doc$nuclear_placed_DOC)
intersect_doc$HetNRHomratio <- intersect_doc$nHets/intersect_doc$nNonRefHom

#Non linear association therefore pearson's correlation isn't useful
#cor.test(intersect_doc$nuclear_placed_DOC,intersect_doc$nvariants, method = "pearson")

#install_github("ProcessMiner/nlcor")
c <- nlcor(intersect_doc$nuclear_placed_DOC,intersect_doc$nvariants, plt=T)
c$cor.estimate # 0.62
c$adjusted.p.value # 0.009
print(c$cor.plot)

####Figure 1
cbPalette <- c("#999999", "#FFCCFF", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6600")
x = ggplot(intersect_doc, aes(x=nuclear_placed_DOC,y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Depth of coverage") + geom_point(aes(color=breed)) + scale_x_continuous(limits = c(0,50))+
  scale_y_continuous(labels=comma, limits = c(0,8000000)) + geom_smooth() +
  scale_color_manual(values=cbPalette) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
        axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"),
        legend.title = element_blank())
save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/November_2020/Figures/Fig1_DOC_nvariants_breed.tiff", x, base_height = 3.5, base_width = 6, dpi = 300)


####GB estimation
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper")

#GB
data = read.table("genetic_burden_by_individual.txt", header=F) # V3 = het, V4 = hom
data$total <- data$V3 + data$V4
data$total <- 0.78*data$total
data$V4 <- 0.78*data$V4
mean(data$total) 
range(data$total)
mean(data$V3)
range(data$V3)
mean(data$V4)

range(data$V4)

mean(data$total[data$V2 == "Arabian"]) 
mean(data$total[data$V2 == "Belgian"])
mean(data$total[data$V2 == "Clydesdale"])
mean(data$total[data$V2 == "Icelandic"])
mean(data$total[data$V2 == "Morgan"])
mean(data$total[data$V2 == "QH"])
mean(data$total[data$V2 == "Shetland"])
mean(data$total[data$V2 == "STB"])
mean(data$total[data$V2 == "TB"])
mean(data$total[data$V2 == "WP"])

#Look for association between breed and genetic burden accounting for DOC
#Add in DOC info
DOC <- read.table("../../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "total")
gb_doc <- merge(data,DOC, by="Sample")
gb_br <- gb_doc %>% 
  filter(!grepl('Other', breed))

fit1 <- (lm(total ~ breed, data=gb_br))
fit2 <- (lm(total ~ breed + nuclear_placed_DOC, data=gb_br))
anova(fit1,fit2)

gb_m <- (lm(total ~ breed + nuclear_placed_DOC,data=gb_br))
n_hom_gb_m <- (lm(hom ~ breed + nuclear_placed_DOC,data=gb_br))
summary(gb_m)
summary(n_hom_gb_m)

#Get EMMEANs
gb_emm <- emmeans(gb_m, "breed", weights = "proportional", type = "response")
gb_emm
n_hom_gb_emm <- emmeans(n_hom_gb_m, "breed", weights = "proportional", type = "response")
n_hom_gb_emm

#GET INFORMATION RE Ne
#Get values from emmeans
gb_emm
n_hom_gb_emm
gb_t <- as.data.frame(table(gb_br$breed))
colnames(gb_t) <- c("breed", "nvariants")
gb_t$nvariants <- c(918,983,927,992,770,"NA",908,871,875,766,899)
gb_t$nvariants <- as.numeric(gb_t$nvariants)
gb_t$nhom <- c(230,245,262,256,172,"NA",203,212,204,165,216)
gb_t$nhom <- as.numeric(gb_t$nhom)

# Add in Ne from Jessica's paper
gb_t$ne1 <- "NA"
gb_t$ne1[gb_t$breed == "Arabian"] <- 346
gb_t$ne1[gb_t$breed == "Belgian"] <- 431
gb_t$ne1[gb_t$breed == "Icelandic"] <- 555
gb_t$ne1[gb_t$breed == "Morgan"] <- 448
gb_t$ne1[gb_t$breed == "QH"] <- 426
gb_t$ne1[gb_t$breed == "STB"] <- 290
gb_t$ne1[gb_t$breed == "TB"] <- 190
gb_t$ne1 <- as.numeric(gb_t$ne1)

cor.test(x=gb_t$ne1,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne1,y=gb_t$nhom, method = "pearson", use='complete.obs')

#Add in effective population size (from Sam's paper)
gb_t$ne <- "NA"
gb_t$ne[gb_t$breed == "Arabian"] <- 3561
gb_t$ne[gb_t$breed == "Belgian"] <- 3570
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "Morgan"] <- 4481
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "QH"] <- 6516
gb_t$ne[gb_t$breed == "STB"] <- 2528
gb_t$ne[gb_t$breed == "TB"] <- 1784
gb_t$ne[gb_t$breed == "WP"] <- 5625

gb_t$ne <- as.numeric(gb_t$ne)
cor.test(x=gb_t$ne,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne,y=gb_t$nhom, method = "pearson", use='complete.obs')


####Figure 2a and b
#Plot EMMEANS
#Number of variants
x <- plot(gb_emm) + geom_boxplot(orientation = 90) + theme_bw() + xlab("EMMEAN of genetic burden") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits=c(700,1050)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#Number of homozygous variants
x1 <- plot(n_hom_gb_emm) + geom_boxplot(orientation = 90) + theme_bw() + xlab("EMMEAN of homozygous genetic burden") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper/")

####LOF variants
data = read.table("lof_by_individual.txt", header=F) # V3 = het, V4 = hom, V5 = missing?
data$total = data$V3 + data$V4
data$total <- 0.78*data$total
data$V4 <- 0.78*data$V4
mean(data$total) 
range(data$total)
mean(data$V3)
range(data$V3)
mean(data$V4)
range(data$V4)

mean(data$total[data$V2 == "Arabian"]) 
mean(data$total[data$V2 == "Belgian"]) 
mean(data$total[data$V2 == "Clydesdale"]) 
mean(data$total[data$V2 == "Icelandic"]) 
mean(data$total[data$V2 == "Morgan"]) 
mean(data$total[data$V2 == "QH"]) 
mean(data$total[data$V2 == "Shetland"]) 
mean(data$total[data$V2 == "STB"]) 
mean(data$total[data$V2 == "TB"]) 
mean(data$total[data$V2 == "WP"]) 


#Look for association between breed and LOF accounting for DOC
#Add in DOC info
DOC <- read.table("../../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "total")
gb_doc <- merge(data,DOC, by="Sample")
gb_br <- gb_doc %>% 
  filter(!grepl('Other', breed))

fit1 <- (lm(total ~ breed, data=gb_br))
fit2 <- (lm(total ~ breed + nuclear_placed_DOC, data=gb_br))
anova(fit1,fit2)

gb_m <- (lm(total ~ breed + nuclear_placed_DOC,data=gb_br))
n_hom_gb_m <- (lm(hom ~ breed + nuclear_placed_DOC,data=gb_br))
summary(gb_m)
summary(n_hom_gb_m)

#Get EMMEANs
gb_emm <- emmeans(gb_m, "breed", weights = "proportional", type = "response")
gb_emm
n_hom_gb_emm <- emmeans(n_hom_gb_m, "breed", weights = "proportional", type = "response")
n_hom_gb_emm

#Get values from emmeans
gb_t <- as.data.frame(table(gb_br$breed))
colnames(gb_t) <- c("breed", "nvariants")
gb_t$nvariants <- c(755,807,752,817,615,"NA",744,706,709,626,729)
gb_t$nvariants <- as.numeric(gb_t$nvariants)
gb_t$nhom <- c(180,192,199,199,127,"NA",159,161,155,128,165)
gb_t$nhom <- as.numeric(gb_t$nhom)

# Add in Ne from Jessica's paper
gb_t$ne1 <- "NA"
gb_t$ne1[gb_t$breed == "Arabian"] <- 346
gb_t$ne1[gb_t$breed == "Belgian"] <- 431
gb_t$ne1[gb_t$breed == "Icelandic"] <- 555
gb_t$ne1[gb_t$breed == "Morgan"] <- 448
gb_t$ne1[gb_t$breed == "QH"] <- 426
gb_t$ne1[gb_t$breed == "STB"] <- 290
gb_t$ne1[gb_t$breed == "TB"] <- 190
gb_t$ne1 <- as.numeric(gb_t$ne1)

cor.test(x=gb_t$ne1,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne1,y=gb_t$nhom, method = "pearson", use='complete.obs')

#Add in effective population size (from Sam's paper)
gb_t$ne <- "NA"
gb_t$ne[gb_t$breed == "Arabian"] <- 3561
gb_t$ne[gb_t$breed == "Belgian"] <- 3570
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "Morgan"] <- 4481
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "QH"] <- 6516
gb_t$ne[gb_t$breed == "STB"] <- 2528
gb_t$ne[gb_t$breed == "TB"] <- 1784
gb_t$ne[gb_t$breed == "WP"] <- 5625

gb_t$ne <- as.numeric(gb_t$ne)
cor.test(x=gb_t$ne,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne,y=gb_t$nhom, method = "pearson", use='complete.obs')

####Figure 2c and d
####Plot EMMEANS
#Number of variants
x3 <- plot(gb_emm) + geom_boxplot(orientation = 90) + theme_bw() + xlab("EMMEAN of LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits = c(540,900)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#Number of homozygous variants
x4 <- plot(n_hom_gb_emm) + geom_boxplot(orientation = 90) + theme_bw() + xlab("EMMEAN of homozygous LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits =c(100,220)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#Combine with GB figures a and b for paper
first_col <- plot_grid(x,x1, labels = c("a", "b"),rel_widths = 1,rel_heights = 1, ncol = 1)
second_col <- plot_grid(x3,x4, labels = c("c", "d"), rel_widths = 1,rel_heights = 1, ncol = 1)
x_c <- plot_grid(first_col,second_col, ncol = 1, rel_widths = 1, rel_heights = 1)

save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/November_2020/Figures/Fig2_gb_gb_hom_EMMEANS.tiff", x_c, base_height = 12, base_width = 6, dpi = 200)

####Known variants
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/")
#Need to restructure my data: disease/breed/genotype/count
data1 <- read.table("variants_by_individual_R.txt", header=T)
deleterious_causative <- data1 %>% 
  filter(!grepl('n', disease)) %>%
  filter(!grepl('n', causative))

bp1 <- ggplot(deleterious_causative, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,8)) + 
  xlab("Causal variants (disease)") +
  scale_fill_manual(values=cbPalette) +
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))

#Print for Havemeyer presentation
bp1 <- ggplot(deleterious_causative, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = "stack") +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,25)) + 
  xlab("Causal variants (disease)") +
  scale_fill_manual(values=cbPalette) +
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))
save_plot("/Users/durwa004/Documents/Postdoc/Projects/Genetic_burden/Havemeyer_2020/DCV_OMIA_variants.tiff", 
          bp1, base_height = 12,base_width = 24, dpi = 100)

non_deleterious_causative <- data1 %>% 
  filter(!grepl('y', disease)) %>%
  filter(!grepl('n', causative))

bp2 <- ggplot(non_deleterious_causative, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual(values = c(0.5,1)) +
  scale_fill_manual(values=cbPalette) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,150)) + 
  xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))

#For Havemeyer
bp2 <- ggplot(non_deleterious_causative, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = "stack") +
  scale_alpha_manual(values = c(0.5,1)) +
  scale_fill_manual(values=cbPalette) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,500)) + 
  xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))
#Save for Havemeyer
save_plot("/Users/durwa004/Documents/Postdoc/Projects/Genetic_burden/Havemeyer_2020/Non_DCV_OMIA_variants.tiff", 
          bp2, base_height = 12,base_width = 24, dpi = 100)

deleterious_assoc <- data1 %>%
  filter(!grepl('y', causative))

bp3 <- ggplot(deleterious_assoc, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,150)) + 
  xlab("Associated variants") + 
  scale_fill_manual(values=cbPalette) +
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))

#For Havemeyer
bp3 <- ggplot(deleterious_assoc, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = "stack") +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,600)) + 
  xlab("Associated variants") + 
  scale_fill_manual(values=cbPalette) +
  theme(axis.text = element_text(size=14,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=16,face="bold"))

save_plot("/Users/durwa004/Documents/Postdoc/Projects/Genetic_burden/Havemeyer_2020/Assoc_OMIA_variants.tiff", 
          bp3, base_height = 12,base_width = 24, dpi = 100)

#Save as combined plot
first_row <- plot_grid(bp1, labels = c("a"),rel_widths = 1)
second_row <- plot_grid(bp2,labels = c("b"), rel_widths = 1)
third_row <- plot_grid(bp3,labels = c("c"), rel_widths = 1)
x_c <- plot_grid(first_row,second_row, third_row,
                 ncol = 1, rel_widths = 1, rel_heights = 1)

#Save as dual plot
save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/November_2020/Figures/Fig3_OMIA_variants.tiff", 
          x_c, base_height = 12,base_width = 24, dpi = 100)


####QTLs
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("QTLs_table.txt", header=T)
data$genotype <- as.factor(data$genotype)

bp1 <- ggplot(data, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Number of horses") + scale_y_continuous() + xlab("QTLs") +
  scale_fill_manual(values=cbPalette) + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"),
        axis.text.x=element_blank())

bp2 <- ggplot(data, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Number of horses") + scale_y_continuous() + xlab("QTLs") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"),
        axis.text.x=element_blank())

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("a", "b"), ncol=1)
#Save as dual plot
save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/November_2020/Figures/SI_1_QTLs.tiff", first_row, base_height = 12,base_width = 24, dpi=100)
