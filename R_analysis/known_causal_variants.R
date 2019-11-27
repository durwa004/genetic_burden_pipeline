library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("known_variants_present.txt", header=T)

bp <- ggplot(data, aes(x=disease, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")


data$Hom_het <- as.factor(data$Hom_het)
bp <- ggplot(data, aes(x=disease, fill=Hom_het)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as dual plot
save_plot("Coat_color_ACs_hom_het.tiff", bp, base_height = 3,base_width = 6)
save_plot("Coat_color_ACs.tiff", bp, base_height = 3,base_width = 6)

data1 = read.table("known_dz_with_breeds.txt", header=T)
bp <- ggplot(data1, aes(x=disease, fill=breed)) + geom_histogram(stat = "count") + 
  ylab("Number of horses") + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold")) + scale_fill_brewer(palette = "Accent")

data1$Hom_het <- as.factor(data1$Hom_het)
bp <- ggplot(data1, aes(x=disease, fill=Hom_het)) +  geom_histogram(stat = "count")  + 
  ylab("Number of horses") + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
save_plot("Dz_variants_ACs.tiff", bp, base_height = 3,base_width = 6)
save_plot("Dz_variants_ACs_hom_het.tiff", bp, base_height = 3,base_width = 6)
