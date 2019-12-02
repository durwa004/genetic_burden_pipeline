library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("non-dz_variants_table.txt", header=T)
data1 = read.table("dz_variants_table.txt", header=T)

bp1 <- ggplot(data, aes(x=Disease, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

data$genotype <- as.factor(data$genotype)
bp2 <- ggplot(data, aes(x=Disease, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")


bp3 <- ggplot(data1, aes(x=Disease, fill=breed)) + geom_histogram(stat = "count") + 
  ylab("Number of horses") + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

data1$genotype <- as.factor(data1$genotype)
bp4 <- ggplot(data1, aes(x=Disease, fill=genotype)) +  geom_histogram(stat = "count")  + 
  ylab("Number of horses") + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"))
second_row <- plot_grid(bp3,bp4, labels = c("C", "D"))
known_variants <- plot_grid(first_row, second_row, ncol =1)
#Save as dual plot
save_plot("known_variants.tiff", known_variants, base_height = 12,base_width = 24)
