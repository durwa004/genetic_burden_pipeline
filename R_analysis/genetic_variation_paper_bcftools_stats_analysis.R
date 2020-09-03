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

#Only include autosomes and chr X (not MT and unplaced contigs)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")

####Figure out number of variants per individual
intersect_stats <- read.table("intersect_by_ind_number_of_variants.txt",header=T)
intersect_stats$nvariants <- intersect_stats$nNonRefHom + intersect_stats$nHets
mean(intersect_stats$nvariants)
mean(intersect_stats$nIndels)
mean(intersect_stats$nNonRefHom)
mean(intersect_stats$nHets)
intersect_stats$tstv <- intersect_stats$Ts/intersect_stats$Tv

#Stuff for Lauren
PPID <- as.data.frame(intersect_stats[175:210,c("Sample", "nvariants")])
write.table(PPID, "/Users/durwa004/Desktop/PPID_nvariants.txt",quote=F)

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

#Add in DOC info
DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) = c("Sample", "total_DOC","nuclear_placed_DOC")
intersect_doc <- merge(intersect_stats,DOC, by="Sample")
summary(intersect_doc$nuclear_placed_DOC)
intersect_doc$HetNRHomratio <- intersect_doc$nHets/intersect_doc$nNonRefHom

table(intersect_doc$breed)

kruskal.test(intersect_doc$HetNRHomratio, intersect_doc$breed)
kruskal.test(intersect_doc$tstv, intersect_doc$nuclear_placed_DOC)

intersect_br <- intersect_doc %>% 
  filter(!grepl('Other', breed))

kruskal.test(intersect_br$HetNRHomratio, intersect_br$breed)
kruskal.test(intersect_br$tstv, intersect_br$breed)

summary(lm(intersect_br$tstv ~ intersect_br$breed))
gb_m <- (lm(tstv ~ breed + nuclear_placed_DOC,data=intersect_br))
n_hom_gb_m <- (lm(HetNRHomratio ~ breed + nuclear_placed_DOC,data=intersect_br))
n_var_m <- (lm(nvariants ~ breed + nuclear_placed_DOC,data=intersect_br))
n_hom_var_m <- (lm(nNonRefHom ~ breed + nuclear_placed_DOC,data=intersect_br))
n_var_m_all <- (lm(nvariants ~ nuclear_placed_DOC,data=intersect_br))
n_hom_var_m_all <- (lm(nNonRefHom ~ nuclear_placed_DOC,data=intersect_br))
n_indel_m_all <- (lm(nIndels ~ nuclear_placed_DOC,data=intersect_br))

summary(gb_m)
summary(n_hom_gb_m)
summary(n_var_m)

#Get EMMEANs
gb_emm <- emmeans(gb_m, "breed", weights = "proportional", type = "response")
gb_emm
test(gb_emm, null = mean(intersect_br$tstv))

n_hom_gb_emm <- emmeans(n_hom_gb_m, "breed", weights = "proportional", type = "response")
n_hom_gb_emm
test(n_hom_gb_emm, null = mean(intersect_br$HetNRHomratio))
n_var_emm <- emmeans(n_var_m, "breed", weights = "proportional", type = "response")
n_var_emm
n_hom_var_emm <- emmeans(n_hom_var_m, "breed", weights = "proportional", type = "response")
n_hom_var_emm
n_var_emm_all <- emmeans(n_var_m_all, "nuclear_placed_DOC", weights = "proportional", type = "response")
n_var_emm_all
n_hom_var_emm_all <- emmeans(n_hom_var_m, "nuclear_placed_DOC", weights = "proportional", type = "response")
n_hom_var_emm_all
n_indel_emm_all <- emmeans(n_indel_m_all, "nuclear_placed_DOC", weights = "proportional", type = "response")
n_indel_emm_all

gb_t <- as.data.frame(table(intersect_br$breed))
colnames(gb_t) <- c("breed", "nvariants")
gb_t$nvariants <- c(5504169,6247468,5834602,6214390,5625180,"NA",5562870,5641106,5612116,4858946,5948175)
gb_t$nvariants <- as.numeric(gb_t$nvariants)
gb_t$nhom <- c(1896660,2178546,2344923,2196847,1804319,"NA",1608857,1875154,1860016,1367011,1950208)
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

#Look for correlation between number of variants and DOC with line of best fit 
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
save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/chapt2_genetic_variation/DOC_nvariants_breed.tiff", x, base_height = 3.5, base_width = 6, dpi = 300)


#split DOC into quartiles
i_stats <- intersect_doc
i_stats$DOC <- with(i_stats,cut(nuclear_placed_DOC, breaks=quantile(nuclear_placed_DOC,
                                                                    probs=seq(0,1,by=0.25)),include.lowest = TRUE))
levels(i_stats$DOC) <- c("Q1", "Q2", "Q3", "Q4")
summary(intersect_doc$nuclear_placed_DOC)

#Show interaction between nvariants and DOC
variants_DOC <- (lm(nvariants ~ breed + DOC,data=i_stats))

x = emmip(variants_DOC, breed ~ DOC, cov.reduce = FALSE, 
          type = "response") + theme_bw() + xlab("DOC quantiles") + 
  ylab("Number of variants") +scale_y_continuous(labels=comma) +  geom_point(aes(color=breed)) +
  scale_color_manual(values=cbPalette) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"),
        legend.title = element_blank()) 
save_plot("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/chapt2_genetic_variation/Relationship_between_nvariants_DOC.tiff", x, base_height = 3.5, base_width = 8)


#Add in breed colors
x = ggplot(intersect_doc, aes(x=nuclear_placed_DOC,y=nvariants)) + theme_bw() + ylab("Number of variants") + 
  xlab("Depth of coverage") + geom_point(aes(color=breed)) + scale_x_continuous(limits = c(0,50))+
  scale_y_continuous(labels=comma, limits = c(0,8000000)) + geom_smooth() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
        axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"),
        legend.title = element_blank())
save_plot("/Users/durwa004/Documents/PhD/Thesis/Thesis/Post_defense_edits/DOC_nvariants_breed.tiff", x, base_height = 3.5, base_width = 6)

fit1 <- (lm(nvariants ~ breed, data=intersect_br))
fit2 <- (lm(nvariants ~ breed + nuclear_placed_DOC, data=intersect_br))
anova(fit1,fit2)

variants_m <- (lm(nvariants ~ breed + nuclear_placed_DOC,data=intersect_br))
variants_h_m <- (lm(nNonRefHom ~ breed + nuclear_placed_DOC,data=intersect_br))

#Get EMMEANs
nvariants_emm <- emmeans(variants_m, specs = "breed", weights = "proportional", 
                         type = "response")
nvariants_hom_emm <- emmeans(variants_h_m, specs = "breed", weights = "proportional", 
                         type = "response")

####Plot EMMEANS
#Number of variants
x1 <- plot(nvariants_emm) + geom_boxplot() + theme_bw() + 
  xlab("EMMEAN of number of variants") + ylab("Breed") + 
  scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
x2 <- plot(nvariants_hom_emm) + geom_boxplot() + theme_bw() + 
  xlab("EMMEAN of number of homozygous variants") + ylab("Breed") + 
  scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#Save as combined plot
first_row <- plot_grid(x1,x2, labels = c("A", "B"), ncol=1)
#Save as dual plot
save_plot("../../Papers_for_publication/Nature_genetics/Figures/nvariants_nhom_EMMEANS.tiff", first_row, base_height = 6,base_width = 12)


#Non linear association therefore pearson's correlation isn't useful
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

#Get number of variants shared by breed figure
#May need to do a Venn diagram
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")
shared = read.table("breed_breed_shared_unique_variants.txt", header=T)

row_means_shared <- data.frame(ID=shared[,1],Means=rowMeans(shared[,2:11], na.rm=T))
mean(row_means_shared$Means)
colMax <- function(data) sapply(data,max,na.rm=TRUE)
colMax(shared[,2:11])
colSort <- function(data, ...) sapply(data, sort, ...)
colSort(shared[2:11], decreasing=TRUE)

library(reshape2)
df$row.names<-rownames(df)
long.df<-melt(df,id=c("row.names"))
plotted<-ggplot(long.df,aes(x=row.names,y=variable,color=value))+geom_point()

x = ggplot(shared, aes(x=breed, y=tstv)) + theme_bw() + ylab("Number of indels") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../bcftools_stats_output/10_breeds_tstv.tiff", x, base_height = 3.5, base_width = 6)



####AF differences by breed - not sure if this works?
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")
#Categorize AF: Using python
AF_data <- read.table("all_horses_AF_freq_info.txt",header=T)
AF_data$means <- rowMeans(AF_data[,2:536])

#AF
x = ggplot(AF_data, aes(x=AF, y=means)) + theme_bw() + ylab("Number of variants") + 
  xlab("Allele frequency (%)") + geom_point() + scale_y_continuous(labels=comma, limits= c(0,3000000)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(angle=90), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
save_plot("535_horses_AF.tiff", x, base_height = 3.5, base_width = 6)

sum(AF_data$means[AF_data$AF < 5]) # 18108154
sum(AF_data$means[AF_data$AF >= 5]) # 13530513
sum(AF_data$means[AF_data$AF < 5]) + sum(AF_data$means[AF_data$AF >= 5]) # 31638667
18108154/31638667
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


#Get EMMEANS of number of variants
#Get regions with more/less variation
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")

#Read in number of variants for each region
regions <- read.table("regions_number_of_variants.txt", header=T)
regions_v <- sum(regions$no_records)
regions_snp <- sum(regions$no_SNPs)
regions_mnp <- sum(regions$no_MNPs)
regions_indel <- sum(regions$no_indels)
regions_ma <- sum(regions$no_multiallelic_sites)
regions_ma_snp <- sum(regions$no_nultiallelic_SNPs)
regions_tstv <- mean(regions$tstv)

mean(regions$no_records)
mean(regions$no_SNPs)
mean(regions$no_indels)
range(regions$no_records)
range(regions$no_SNPs)
range(regions$no_indels)
#High variation regions:
regions_high <- regions %>% filter(no_records > (mean(regions$no_records))*2)
str(regions_high)
mean(regions_high$no_records)
mean(regions_high$no_SNPs)
mean(regions_high$no_indels)
mean(regions_high$tstv)
#Print out high variation chrom/pos1/pos2 regions
write.table(regions_high, file = "High_variation_regions.txt", quote=F,sep = "\t")

#Low variation regions:
regions_low <- regions %>% filter(no_records < (mean(regions$no_records))/2)
str(regions_low)
mean(regions_low$no_records)
mean(regions_low$no_SNPs)
mean(regions_low$no_indels)
mean(regions_low$tstv)
#Print out low variation chrom/pos1/pos2 regions
write.table(regions_low, file = "Low_variation_regions.txt", quote=F,sep = "\t")


