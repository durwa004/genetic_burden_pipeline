library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("known_variants_with_breed.txt", header=T)
data$genotype <- as.factor(data$genotype)

deleterious_causative <- data %>% 
  filter(!grepl('n', deleterious)) %>%
  filter(!grepl('n', causative))

bp1 <- ggplot(deleterious_causative, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits= c(0,25)) + xlab("Causal variants (deleterious)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

bp2 <- ggplot(deleterious_causative, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits= c(0,25)) + xlab("Causal variants (deleterious)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"))
#Save as dual plot
save_plot("causative_deleterious.tiff", first_row, base_height = 12,base_width = 24)


non_deleterious_causative <- data %>% 
  filter(!grepl('y', deleterious)) %>%
  filter(!grepl('n', causative))

bp3 <- ggplot(non_deleterious_causative, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count") + 
  ylab("Number of horses") + xlab("Causal variants (non-deleterious)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

bp4 <- ggplot(non_deleterious_causative, aes(x=Phenotype, fill=genotype)) +  geom_histogram(stat = "count")  + 
  ylab("Number of horses") + xlab("Causal variants (non-deleterious)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
second_row <- plot_grid(bp3,bp4, labels = c("A", "B"))
#Save as dual plot
save_plot("causative_non_deleterious.tiff", second_row, base_height = 12,base_width = 24)

assoc_variants <- data %>%
  filter(!grepl('y', causative))
bp1 <- ggplot(deleterious_assoc, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous() + xlab("Associated variants") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

bp2 <- ggplot(deleterious_assoc, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous() + xlab("Associated variants") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"))
#Save as dual plot
save_plot("associated_variants.tiff", first_row, base_height = 12,base_width = 24)


