library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

setwd("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/")
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
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"))

non_deleterious_causative <- data1 %>% 
  filter(!grepl('y', disease)) %>%
  filter(!grepl('n', causative))

bp2 <- ggplot(non_deleterious_causative, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual(values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,150)) + 
  xlab("Causal variants (non-disease)") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"))

deleterious_assoc <- data1 %>%
  filter(!grepl('y', causative))

bp3 <- ggplot(deleterious_assoc, aes(x = Phenotype, y = genotype_count, fill = breed, group = genotype)) + 
  geom_bar(stat = "identity", aes(alpha = factor(genotype)),
           position = position_dodge(width = 1)) +
  scale_alpha_manual("Genotype", values = c(0.5,1)) +
  ylab("Genotype Count") + scale_y_continuous(limits = c(0,150)) + 
  xlab("Associated variants") + 
  theme(axis.text = element_text(size=10,angle=90,hjust=1),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size=12,face="bold"))

#Save as combined plot
first_row <- plot_grid(bp1, labels = c("A"),rel_widths = 1)
second_row <- plot_grid(bp2,labels = c("B"), rel_widths = 1)
third_row <- plot_grid(bp3,labels = c("C"), rel_widths = 1)
x_c <- plot_grid(first_row,second_row, third_row,
                 ncol = 1, rel_widths = 1, rel_heights = 1)

#Save as dual plot
save_plot("../Draft_May_14/OMIA_variants.tiff", 
          x_c, base_height = 12,base_width = 24)

#Get stats for paper
data_info <- read.table("known_variant_all_tidy.txt", header = T)
deleterious_causative <- data_info %>% 
  filter(!grepl('n', Deleterious)) %>%
  filter(!grepl('n', Causative))

mean(deleterious_causative$AF)
range(deleterious_causative$AF)
length(deleterious_causative$AF)
table(deleterious_causative$AC)

deleterious_associated <- data_info %>% 
  filter(!grepl('n', Deleterious)) %>%
  filter(!grepl('y', Causative))

mean(deleterious_associated$AF)
range(deleterious_associated$AF)
length(deleterious_associated$AF)
table(deleterious_associated$AC)

non_deleterious_causative <- data_info %>% 
  filter(!grepl('y', Deleterious)) %>%
  filter(!grepl('n', Causative))

mean(non_deleterious_causative$AF)
range(non_deleterious_causative$AF)
length(non_deleterious_causative$AF)
table(non_deleterious_causative$AC)

non_deleterious_associated <- data_info %>% 
  filter(!grepl('y', Deleterious)) %>%
  filter(!grepl('y', Causative))

mean(non_deleterious_associated$AF)
range(non_deleterious_associated$AF)
length(non_deleterious_associated$AF)
table(non_deleterious_associated$AC)

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

################################################################################
###############################################################################
#QTLs
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/")
data = read.table("QTLs_table.txt", header=T)
data$genotype <- as.factor(data$genotype)

bp1 <- ggplot(data, aes(x=Phenotype, fill=breed)) + geom_histogram(stat = "count")  + 
  ylab("Number of horses") + scale_y_continuous() + xlab("QTLs") + 
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
first_row <- plot_grid(bp1,bp2, labels = c("A", "B"), ncol=1)
#Save as dual plot
save_plot("../../Papers_for_publication/Nature_genetics/Figures/QTLs.tiff", first_row, base_height = 12,base_width = 24)

