library(scales)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(reshape2)
library(dplyr)

test <- read.table("/Users/durwa004/Downloads/genetic_burden_details.txt", header = T)
test2 <- test[,1:2]
test$pos <- NULL
test$chrom <- NULL
ID_test <- colnames(test)
het_test <- data.frame(colSums(test == "het"))
colnames(het_test) <- "het_test"
hom_test <- as.data.frame(colSums(test == "hom"))
colnames(hom_test) <- "hom_test"

test.df <- cbind(ID_test, het_test$het_test, hom_test$hom_test)
test.df <- as.data.frame(test.df)
colnames(test.df) <- c("ID", "het_test", "hom_test")
test.df$het_test <- as.numeric(test.df$het_test)
test.df$hom_test <- as.numeric(test.df$hom_test)

test.df$GB <- test.df$het_test + test.df$hom_test

####Look for outliers
test.df$hom_l <- log10(test.df$hom_test)
test.df$het_l <- log10(test.df$het_test)

lower_bound <- quantile(test.df$het_test, 0.01)
upper_bound <- quantile(test.df$het_test, 0.99)
highlight_df1 <- test.df %>%
  filter(het_test > upper_bound)
highlight_df2 <- test.df %>%
  filter(het_test < lower_bound)

test.df %>% 
  ggplot(aes(x=ID, y = het_test))+ 
  scale_color_brewer(palette = "Paired") + 
  geom_point(alpha=0.3) + 
  geom_point(data=highlight_df1, 
             aes(x=ID, y = het_test, 
                 color= ID),
             size=3)+ 
  geom_point(data=highlight_df2, 
             aes(x=ID, y = het_test, 
                 color= ID),
             size=3) 


lower_bound <- quantile(test.df$hom_test, 0.01)
upper_bound <- quantile(test.df$hom_test, 0.99)
highlight_df1 <- test.df %>%
  filter(hom_test > upper_bound)
highlight_df2 <- test.df %>%
  filter(hom_test < lower_bound)

test.df %>% 
  ggplot(aes(x=ID, y = hom_test)) + 
  geom_point(alpha=0.3) + 
  geom_point(data=highlight_df1, 
             aes(x=ID, y = hom_test, 
                 color= ID),
             size=3)







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

