library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)
library(ggrepel)
library(Gviz)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/")
lof <- read.table("lof_details_brief.txt", header=T,sep="\t")

#Have a look at genes
y <- as.data.frame(table(lof$gene_ann))
colnames(y) <- c("Gene", "count")

#Get labels for genes with over 10 LOF variants
y %>% filter(y$count > 10)
y$fill <- ifelse((y$count > 10), "black", "red")
#Histogram of number of LOF variants per gene
x = ggplot(y, aes(x = Gene, y = count, label = Gene, fill = fill)) +
  geom_histogram(stat = "identity") + ylab("Number of LOF variants") + 
  xlab("Genes containing LOF variants") + scale_y_continuous(limits = c(0,30)) + 
  geom_text_repel(data = subset(y, count > 10)) + scale_fill_manual(values = c("red", "black")) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_per_gene.jpeg", x, base_height = 3, base_width = 6)

#Histogram of number of genes
y_s <- y[(y$count > 5),]
y_s$Gene <- as.character(y_s$Gene)
y_s$Gene <- as.factor(y_s$Gene)

x = ggplot(y_s, aes(x = count)) +
  geom_histogram(binwidth=1, color="pink",fill="grey") + ylab("Frequency") + 
  xlab("Number of LOF\nvariants per gene") + scale_y_continuous(limits = c(0,20)) + 
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_multiple_LOF_per_gene.jpeg", x, base_height = 3, base_width = 6)

#Do for high AF genes
se_e <- droplevels(lof %>% group_by(gene_ann) %>% summarise(se_e = mean(AF,na.rm=T)))
colnames(se_e) <- c("gene", "AF")
se_e$fill <- ifelse((se_e$AF > 0.98), "black", "red")
se_e <- as.data.frame(se_e)
x = ggplot(se_e, aes(x = gene, y = AF, label = gene, fill = fill)) +
  geom_bar(stat = "identity") + ylab("Allele frequency") + 
  xlab("Genes containing LOF variants") + scale_y_continuous(limits = c(0,1)) + 
  geom_text_repel(data = subset(se_e, AF > 0.98)) + scale_fill_manual(values = c("red", "black")) +
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_AF_per_gene.jpeg", x, base_height = 3, base_width = 6)

high_AF <- droplevels(lof %>% filter(AF > 0.98))
y %>% filter(y$count > 10)
y$fill <- ifelse((y$count > 10), "black", "red")

high_AF$gene_ann <- as.character(high_AF$gene_ann)
high_AF$gene_ann <- as.factor(high_AF$gene_ann)



