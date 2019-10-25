library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)
library(scales)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("known_CC_with_breeds.txt", header=T)

bp <- ggplot(data, aes(x=disease, y=AC, fill=breed)) + geom_bar(stat="identity") + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), axis.title = element_text(size=12,face="bold"))
save_plot("Coat_color_ACs.tiff", bp, base_height = 3,base_width = 6)

df2 <- as.data.frame(data$disease)
df2$AC <- data$AC
df2$Hom_het <- data$Hom_het
df3 <- melt(df2, id=c("data$disease","AC"))
df4 <- df3[order(df3$value),]
bp <- ggplot(df4, aes(x=data$disease, y=AC, fill=value)) + geom_bar(stat="identity") + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), axis.title = element_text(size=12,face="bold"))
save_plot("Coat_color_ACs_hom_het.tiff", bp, base_height = 3,base_width = 6)

data1 = read.table("known_dz_with_breeds.txt", header=T)
bp <- ggplot(data1, aes(x=disease, y=AC, fill=breed)) + geom_bar(stat="identity") + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), axis.title = element_text(size=12,face="bold"))
save_plot("Dz_variants_ACs.tiff", bp, base_height = 3,base_width = 6)

df2 <- as.data.frame(data1$disease)
df2$AC <- data1$AC
df2$Hom_het <- data1$Hom_het
df3 <- melt(df2, id=c("data1$disease","AC"))
df4 <- df3[order(df3$value),]
bp <- ggplot(df4, aes(x=data1$disease, y=AC, fill=value)) + geom_bar(stat="identity") + 
  ylab("Allele Count") + scale_y_continuous(labels=comma) + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), axis.title = element_text(size=12,face="bold"))
save_plot("Dz_variants_ACs_hom_het.tiff", bp, base_height = 3,base_width = 6)


