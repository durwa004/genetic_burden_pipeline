library(ggplot2)
library(cowplot)
library(dplyr)

#Only include autosomes and chr X (not MT and unplaced contigs)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/R_analysis/")
###bcftools
bcftools <- read.table("bcftools_number_of_variants.txt", header=T)
bcf_snp <- sum(bcftools$no_SNPs)
bcf_indel <- sum(bcftools$no_indels)
bcf_tstv <- mean(bcftools$tstv)

###gatk
gatk <- read.table("gatk_number_of_variants.txt", header=T)
gatk_snp <- sum(gatk$no_SNPs)
gatk_indel <- sum(gatk$no_indels)
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

#Create venn diagram
source("http://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph"))
library(devtools)
install_github("js229/Vennerable")
library(Vennerable)

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

#Number of variants by chromosome length
x = ggplot(intersect, aes(x=CHROM, y=variant_ratio)) + theme_bw() + ylab("Variant:chr length") + 
  xlab("Chromosome") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("intersect_variants.tiff", x, base_height = 3.5, base_width = 6)

####Figure out number of variants per individual
intersect_stats <- read.table("../../bcftools_stats_output/intersect_by_ind_number_of_variants.txt",header=T)
intersect_stats$nvariants <- intersect_stats$nNonRefHom + intersect_stats$nHets
mean(intersect_stats$nvariants)
intersect_stats$tstv <- intersect_stats$Ts/intersect_stats$Tv

####Figure out differences in the number of variants per breed
intersect_stats_br <- intersect_stats %>% 
  filter(!grepl('Other', breed))
kruskal.test(intersect_stats_br$nvariants, intersect_stats_br$breed)
x = ggplot(intersect_stats_br, aes(x=breed, y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("10_breeds_nvariants.tiff", x, base_height = 3.5, base_width = 6)

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
save_plot("10_breeds_nindels.tiff", x, base_height = 3.5, base_width = 6)

###TsTv
kruskal.test(intersect_stats_br$tstv, intersect_stats_br$breed)
mean(intersect_stats$tstv)
x = ggplot(intersect_stats_br, aes(x=breed, y=tstv)) + theme_bw() + ylab("Number of indels") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("10_breeds_tstv.tiff", x, base_height = 3.5, base_width = 6)


####AF differences by breed
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")
#Categorize AF: Using python
AF_data <- read.table("all_horses_AF_freq_info.txt",header=T)
AF_data$means <- rowMeans(AF_data)

af <- read.table("AF_freq_files/M1005_STB_AF_freq.txt", header = T)
AF_data$AF <- af$AF

#AF
x = ggplot(AF_data, aes(x=AF, y=means)) + theme_bw() + ylab("Number of variants") + 
  xlab("AF category") + geom_point(fill=rgb(122/255,0/255,25/255,1)) + scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("548_horses_AF.tiff", x, base_height = 3.5, base_width = 6)

#AF categories
categories <- read.table("all_horses_AF_categorized.txt", header= T)
categories$means <- AF_data$means
x = ggplot(categories, aes(x=AF, y=means)) + theme_bw() + ylab("Number of variants") + 
  xlab("AF category") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) + scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("548_horses_AF_categories.tiff", x, base_height = 3.5, base_width = 6)


#Mean number of singletons
singleton <- AF_data[1,2:548]
singleton <- as.numeric(singleton)
mean(singleton[1:547])

singleton_data <- read.table("10_breeds_singletons.txt",header=F)
kruskal.test(singleton_data$V3, singleton_data$V2)
x = ggplot(singleton_data, aes(x=V2, y=V3)) + theme_bw() + ylab("Number of singletons") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("10_breeds_singletons.tiff", x, base_height = 3.5, base_width = 6)
