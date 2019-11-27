library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)
library(devtools)
library(ggpubr)
library(emmeans)

#Only include autosomes and chr X (not MT and unplaced contigs)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output//")
###bcftools
bcftools <- read.table("bcftools_number_of_variants.txt", header=T)
bcf_v <- sum(bcftools$no_records)
bcf_snp <- sum(bcftools$no_SNPs)
bcf_mnp <- sum(bcftools$no_MNPs)
bcf_indel <- sum(bcftools$no_indels)
bcf_ma <- sum(bcftools$no_multiallelic_sites)
bcf_ma_snp <- sum(bcftools$no_nultiallelic_SNPs)
bcf_tstv <- mean(bcftools$tstv)

###gatk
gatk <- read.table("gatk_number_of_variants.txt", header=T)
gatk_v <- sum(gatk$no_records)
gatk_snp <- sum(gatk$no_SNPs)
gatk_mnp <- sum(gatk$no_MNPs)
gatk_indel <- sum(gatk$no_indels)
gatk_ma <- sum(gatk$no_multiallelic_sites)
gatk_ma_snp <- sum(gatk$no_nultiallelic_SNPs)
gatk_tstv <- mean(gatk$tstv)

###union
union <- read.table("union_number_of_variants.txt", header=T)
union_snp <- sum(union$no_SNPs)
union_indel <- sum(union$no_indels)
union_tstv <- mean(union$tstv)

###intersect
intersect <- read.table("intersect_number_of_variants.txt", header=T)
sum(intersect$no_records)
intersect_snp <- sum(intersect$no_SNPs)
intersect_indel <- sum(intersect$no_indels)
intersect_tstv <- mean(intersect$tstv)
intersect$variant_ratio <- intersect$no_records/intersect$chrom_length
intersect$snp_ratio <- intersect$no_SNPs/intersect$chrom_length
intersect$indel_ratio <- intersect$no_indels/intersect$chrom_length

#Number of variants by chromosome length
x = ggplot(intersect, aes(x=CHROM, y=variant_ratio)) + theme_bw() + ylab("Variant:chr length") + 
  xlab("Chromosome") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) + 
  scale_x_discrete(labels = c(1:31, "X")) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("intersect_variants.tiff", x, base_height = 3.5, base_width = 6)
x = ggplot(intersect, aes(x=CHROM, y=snp_ratio)) + theme_bw() + ylab("SNP:chr length") + 
  xlab("Chromosome") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("intersect_SNPs.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(intersect, aes(x=CHROM, y=indel_ratio)) + theme_bw() + ylab("indel:chr length") + 
  xlab("Chromosome") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("intersect_indels.tiff", x, base_height = 3.5, base_width = 6)

#Create venn diagram
source("http://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph"))
library(devtools)
install_github("js229/Vennerable")
library(Vennerable)

#Get plot for genetic burden
vcombo <- Venn(SetNames=c("ANNOVAR", "SnpEff"),Weight = c(0,2999,735,4993))
jpeg("ANNOVAR_SnpEff_venn.jpeg",width=6,height=6,units="in",res=1350)
plot(vcombo,show=list(SetLabels=FALSE,FaceText=FALSE, Faces=FALSE))
dev.off()

#SNPs
vcombo <- Venn(SetNames=c("gatk","bcftools"),Weight=c(0,gatk_snp-intersect_snp,bcf_snp-intersect_snp,intersect_snp))
plot(vcombo)

##This works to remove text:

jpeg("HC_bcftools_intersect_venn.jpeg",width=6,height=6,units="in",res=1350)
plot(vcombo,show=list(SetLabels=FALSE,FaceText=FALSE,Faces=FALSE))
dev.off()

#indels
vcombo <- Venn(SetNames=c("gatk","bcftools"),Weight=c(0,gatk_indel-intersect_indel,bcf_indel-intersect_indel,intersect_indel))
plot(vcombo)

##This works to remove text:
jpeg("HC_bcftools_intersect_indel_venn.jpeg",width=6,height=6,units="in",res=1350)
plot(vcombo,show=list(SetLabels=FALSE,FaceText=FALSE,Faces=FALSE))
dev.off()


####Figure out number of variants per individual
intersect_stats <- read.table("../with_prze/bcftools_stats_output_with_prze/intersect_by_ind_number_of_variants.txt",header=T)
intersect_stats$nvariants <- intersect_stats$nNonRefHom + intersect_stats$nHets
mean(intersect_stats$nvariants)
intersect_stats$tstv <- intersect_stats$Ts/intersect_stats$Tv

#Add in DOC info
DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) = c("Sample", "total_DOC","nuclear_placed_DOC")
intersect_doc <- merge(intersect_stats,DOC, by="Sample")
summary(intersect_doc$nuclear_placed_DOC)
intersect_doc$HetNRHomratio <- intersect_doc$nHets/intersect_doc$nNonRefHom

table(intersect_doc$breed)
#Details by breed
sum(intersect_doc[intersect_doc$breed == "Arabian",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Arabian",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Arabian",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Arabian",]$tstv)

mean(intersect_doc[intersect_doc$breed == "Belgian",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Belgian",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Belgian",]$tstv)

mean(intersect_doc[intersect_doc$breed == "Clydesdale",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Clydesdale",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Clydesdale",]$tstv)

mean(intersect_doc[intersect_doc$breed == "Icelandic",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Icelandic",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Icelandic",]$tstv)

mean(intersect_doc[intersect_doc$breed == "Morgan",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Morgan",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Morgan",]$tstv)

mean(intersect_doc[intersect_doc$breed == "QH",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "QH",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "QH",]$tstv)

mean(intersect_doc[intersect_doc$breed == "Shetland",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "Shetland",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "Shetland",]$tstv)

mean(intersect_doc[intersect_doc$breed == "STB",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "STB",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "STB",]$tstv)

mean(intersect_doc[intersect_doc$breed == "TB",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "TB",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "TB",]$tstv)

mean(intersect_doc[intersect_doc$breed == "WP",]$nvariants)
mean(intersect_doc[intersect_doc$breed == "WP",]$HetNRHomratio)
mean(intersect_doc[intersect_doc$breed == "WP",]$tstv)

kruskal.test(intersect_doc$HetNRHomratio, intersect_doc$breed)

#Plot DOC histogram
x = ggplot(intersect_doc, aes(x=nuclear_placed_DOC)) + theme_bw() + ylab("Frequency") + 
  xlab("Depth of coverage") + geom_histogram() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/535_individuals_DOC.tiff", x, base_height = 3.5, base_width = 6)

#Look for correlation between number of variants and DOC with line of best fit 
x = ggplot(intersect_doc, aes(x=nuclear_placed_DOC,y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Depth of coverage") + geom_point(fill=rgb(122/255,0/255,25/255,1)) +
  scale_y_continuous(labels=comma, limits = c(0,8000000)) + geom_smooth() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/535_individuals_DOC_nvariants.tiff", x, base_height = 3.5, base_width = 6)

#Non linear association therefore pearson's correlation isn't useful
library(tidyverse)
library(caret)
#Figure out how many orders to use (change the number to max that are signifianct)
model <- lm(nvariants ~ poly(nuclear_placed_DOC,5,raw=TRUE), data = intersect_doc)
predictions <- model %/% predict(test.data)

#cor.test(intersect_doc$nuclear_placed_DOC,intersect_doc$nvariants, method = "pearson")
library(devtools)
#install_github("ProcessMiner/nlcor")
library(nlcor)
c <- nlcor(intersect_doc$nuclear_placed_DOC,intersect_doc$nvariants, plt=T)
c$cor.estimate # 0.62
c$adjusted.p.value # 0.009
print(c$cor.plot)

#Calculate threshold for cut off (cost/benefit type analysis for DOC vs nvariants)

xy = ggplot(intersect_doc, aes(x=nuclear_placed_DOC,y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Depth of coverage") + geom_point(fill=rgb(122/255,0/255,25/255,1)) + geom_smooth()
  scale_y_continuous(labels=comma, limits = c(0,8000000)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

####Figure out differences in the number of variants per breed
intersect_stats_br <- intersect_doc %>% 
  filter(!grepl('Other', breed))
fit1 <- (lm(nvariants ~ breed, data=intersect_stats_br))
fit2 <- (lm(nvariants ~ breed + log10(nuclear_placed_DOC), data=intersect_stats_br))
anova(fit1,fit2)

nvariants_m <- (lm(nvariants ~ breed + log10(nuclear_placed_DOC),data=intersect_stats_br))
summary(nvariants_m)
test(emmeans(nvariants_m, ~breed, "breed",weights = "proportional", type = "response"))
confint(emmeans(nvariants_m, ~breed, "breed",weights = "proportional", type = "response"))

emmip(noise.lm, type ~ size | side)
nvariants_lm_plot <- plotemms(nvariants_m, digits=2, color = "Fat_g") + ylab("Estimated Marginal Mean of log10 TEQ")


kruskal.test(intersect_stats_br$nvariants, intersect_stats_br$breed)
x = ggplot(intersect_stats_br, aes(x=breed, y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/10_breeds_nvariants.tiff", x, base_height = 3.5, base_width = 6)

#DOC
kruskal.test(intersect_stats_br$nuclear_placed_DOC, intersect_stats_br$breed)
x = ggplot(intersect_stats_br, aes(x=breed, y=nuclear_placed_DOC)) + theme_bw() + ylab("Depth of coverage") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/10_breeds_DOC.tiff", x, base_height = 3.5, base_width = 6)

#### Number of indels per individual
mean(intersect_stats$nIndels)

####Figure out differences in the number of variants per breed
intersect_stats_br <- intersect_stats %>% 
  filter(!grepl('Other', breed))
kruskal.test(intersect_stats_br$nIndels, intersect_stats_br$breed)
x = ggplot(intersect_stats_br, aes(x=breed, y=nIndels)) + theme_bw() + ylab("Number of indels") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/10_breeds_nindels.tiff", x, base_height = 3.5, base_width = 6)

###TsTv
kruskal.test(intersect_stats_br$tstv, intersect_stats_br$breed)
mean(intersect_stats$tstv)
x = ggplot(intersect_stats_br, aes(x=breed, y=tstv)) + theme_bw() + ylab("Number of indels") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/10_breeds_tstv.tiff", x, base_height = 3.5, base_width = 6)



####AF differences by breed
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")
#Categorize AF: Using python
AF_data <- read.table("all_horses_AF_freq_info.txt",header=T)
AF_data$means <- rowMeans(AF_data)

af <- read.table("AF_freq_files/M1005_STB_AF_freq.txt", header = T)
AF_data$AF <- af$AF

#AF
x = ggplot(AF_data, aes(x=AF, y=means)) + theme_bw() + ylab("Number of variants") + 
  xlab("Allele frequency") + geom_point(fill=rgb(122/255,0/255,25/255,1)) + scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("535_horses_AF.tiff", x, base_height = 3.5, base_width = 6)

#AF categories
categories <- read.table("all_horses_AF_categorized.txt", header= T)
categories$means <- AF_data$means
categories$AF <- factor(categories$AF,levels = c("<1%", "1-<5%", "5-<10%", "10-<25%", "25-<50%","50-<75%", "75-<100%"))
x = ggplot(categories, aes(x = AF, y = means)) + theme_bw() + ylab("Number of variants") + 
  xlab("AF category") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) + scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("535_horses_AF_categories.tiff", x, base_height = 3.5, base_width = 6)


#Mean number of singletons
singleton <- AF_data[1,2:53]
singleton <- as.numeric(singleton)
mean(singleton)

singleton_data <- read.table("10_breeds_singletons.txt",header=F)
kruskal.test(singleton_data$V3, singleton_data$V2)
x = ggplot(singleton_data, aes(x=V2, y=V3)) + theme_bw() + ylab("Number of singletons") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("10_breeds_singletons.tiff", x, base_height = 3.5, base_width = 6)

#Mean number of reads
read <- read.table("../DOC/reads_by_horse.txt")
#Mean read length
summary(read$V2)
#Mean number of pairs that mapped
summary(read$V3)
