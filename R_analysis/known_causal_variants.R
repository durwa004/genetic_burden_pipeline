library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/")
data = read.table("variants_bt_indvidual.txt", header=T)
data$genotype <- as.factor(data$genotype)

deleterious_causative <- data %>% 
  filter(!grepl('n', disease)) %>%
  filter(!grepl('n', causative))

bp1 <- ggplot(deleterious_causative, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits = c(0,25)) + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

bp2 <- ggplot(deleterious_causative, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits= c(0,25)) + xlab("Causal variants (disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

non_deleterious_causative <- data %>% 
  filter(!grepl('y', disease)) %>%
  filter(!grepl('n', causative))

bp3 <- ggplot(non_deleterious_causative, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count") + 
  ylab("Number of horses") + scale_y_continuous(limits=c(0,500)) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"))

bp4 <- ggplot(non_deleterious_causative, aes(x=Phenotype, fill=genotype)) +  geom_histogram(stat = "count")  + 
  ylab("Number of horses") + scale_y_continuous(limits=c(0,500)) + xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1), 
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"),rel_widths = 1)
second_row <- plot_grid(bp3,bp4, labels = c("C", "D"), rel_widths = 1)
x_c <- plot_grid(first_row,second_row, ncol = 1, rel_widths = 1, rel_heights = 1)

#Save as dual plot
save_plot("../../Paper_2019/Nature_genetics/Useful_figures/causative_variants.tiff", 
          x_c, base_height = 12,base_width = 24)

deleterious_assoc <- data %>%
  filter(!grepl('y', causative))

bp1 <- ggplot(deleterious_assoc, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits = c(0,600)) + xlab("Associated variants") + 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1),
        axis.text.y = element_text(size=10,hjust=1),
        panel.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "in"), legend.key.width = unit(0.3,"in"),
        axis.title = element_text(size=10,face="bold"))

bp2 <- ggplot(deleterious_assoc, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous(limits = c(0,600)) + xlab("Associated variants") + 
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1),
        axis.text.y = element_text(size=10,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"))+ scale_fill_brewer(palette = "Accent")

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"), ncol = 1)
#Save as dual plot
save_plot("../../Paper_2019/Nature_genetics/Useful_figures/associated_variants.tiff", first_row, base_height = 12,base_width = 24)

################################################################################
###############################################################################
#QTLs
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("QTLs_table.txt", header=T)
data$genotype <- as.factor(data$genotype)

bp1 <- ggplot(data, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous() + xlab("QTLs") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"),
        axis.text.x=element_blank())

bp2 <- ggplot(data, aes(x=Phenotype, fill=genotype)) + geom_histogram(stat = "count")  + 
  ylab("Allele Count") + scale_y_continuous() + xlab("QTLs") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.position="none",
        axis.title = element_text(size=12,face="bold"),
        axis.text.x=element_blank())

#Save as combined plot
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"), ncol=1)
#Save as dual plot
save_plot("../Paper_2019/Nature_genetics/Useful_figures/QTLs.tiff", first_row, base_height = 12,base_width = 24)
