library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)
library(ggrepel)
library(Gviz)
library(stringr)

#setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper/")
setwd("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/")
gb_genes<- read.table("genetic_burden_genes_details.txt", header=T,sep="\t")
mean(gb_genes$mean_AF)
range(gb_genes$mean_AF)
mean(gb_genes$n_variants)
range(gb_genes$n_variants)
table(gb_genes$n_variants)
length(gb_genes$gene[gb_genes$n_variants >5])
length(gb_genes$gene[gb_genes$mean_AF >0.5])

m_variants <- (gb_genes$gene[gb_genes$n_variants >5])
sum(str_count(m_variants, "LOC"))

c_variants <- (gb_genes$gene[gb_genes$mean_AF >0.5])
sum(str_count(c_variants, "LOC"))

#Get labels for genes with over 10 LOF variants
y = gb_genes %>% filter(gb_genes$n_variants > 10)
gb_genes$fill <- ifelse((gb_genes$n_variants > 10), "black", "red")

#Histogram of number of LOF variants per gene
x = ggplot(gb_genes, aes(x = gene, y = n_variants, label = gene, fill = fill)) +
  geom_histogram(stat = "identity") + ylab("Number of LOF variants") + 
  xlab("Genes containing LOF variants") + scale_y_continuous(limits = c(0,30)) + 
  geom_text_repel(data = subset(y, n_variants > 10)) + scale_fill_manual(values = c("red", "black")) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
#save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_per_gene.jpeg", x, base_height = 3, base_width = 6)

#Do for high AF genes
xy = gb_genes %>% filter(gb_genes$mean_AF > 0.5)
gb_genes$fill_AF <- ifelse((gb_genes$mean_AF > 0.50), "black", "red")

x1 = ggplot(gb_genes, aes(x = gene, y = mean_AF, label = gene, fill = fill_AF)) +
  geom_bar(stat = "identity") + ylab("Allele frequency") + 
  xlab("Genes containing LOF variants") + scale_y_continuous(limits = c(0,1)) + 
  geom_text_repel(data = subset(gb_genes, mean_AF > 0.95)) + scale_fill_manual(values = c("red", "black")) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
#save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_AF_per_gene.jpeg", x, base_height = 3, base_width = 6)


