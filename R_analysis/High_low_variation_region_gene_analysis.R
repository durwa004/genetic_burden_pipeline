library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)
library(ggrepel)
library(Gviz)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/")
high_genes<- read.table("High_variation_regions_all_brief.txt", header=T,sep="\t")
mean(high_genes$AF)
range(high_genes$AF)
table(high_genes$impact)
table(high_genes$SNP)
table(high_genes$consequence)
high_genes$group = "HIGH"
high_genes$group <- as.factor(high_genes$group)

low_genes<- read.table("Low_variation_regions_all_brief.txt", header=T,sep="\t")
mean(low_genes$AF)
range(low_genes$AF)
table(low_genes$impact)
table(low_genes$SNP)
table(low_genes$consequence)


#Eva recommended converting the total numbers to % to make it easier to see
low_c <- as.data.frame(table(low_genes$consequence))
low_c$percent <- (low_c$Freq/(length(low_genes$consequence)))*100
low_c$group = "LOW"
low_c$group <- as.factor(low_c$group)

high_c <- as.data.frame(table(high_genes$consequence))
high_c$percent <- (high_c$Freq/(length(high_genes$consequence)))*100
high_c$group = "HIGH"
high_c$group <- as.factor(high_c$group)

se_ann <- rbind(low_c, high_c)

consequence_bar <- ggplot(se_ann, aes(x = Var1, y=percent, fill = group,color = group)) +
  geom_histogram(stat = "identity",  position = "dodge2", width = 0.8, alpha = 1) + 
  ylab("Percentage of\ngenic variants") + 
  xlab("Type of variant") +  scale_y_continuous(labels=comma) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key.size = unit(0.1, "in"), legend.key.width = unit(0.1,"in"),
        legend.position = "bottom",
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90))

save_plot("/Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/chapt2_genetic_variation/Figures/high_low_consequence.jpeg", consequence_bar,
          base_height = 4, base_width = 8)

write.table(high_genes$gene, file = "High_variation_regions_gene_symbols.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(low_genes$gene, file = "Low_variation_regions_gene_symbols.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

#Compare AF
#Compare the means
t.test(high_genes$AF, low_genes$AF)

#Just coding regions
high_constraint<- read.table("High_variation_regions_regions_brief_constraint.txt", header=F,sep="\t")
low_constraint<- read.table("Low_variation_regions_regions_brief_constraint.txt", header=F,sep="\t")

mean(high_constraint$V2, na.rm=T)
mean(low_constraint$V2, na.rm=T)
mean(high_constraint$V3, na.rm=T)
mean(low_constraint$V3, na.rm=T)

sum(is.na(high_constraint$V2)) #1296
sum(is.na(high_constraint$V3)) #1305
sum(is.na(low_constraint$V2)) #574
sum(is.na(low_constraint$V3)) #600
length(low_constraint$V1)
length(high_constraint$V1)

t.test(is.na(high_constraint$V2), is.na(low_constraint$V2))
t.test(high_constraint$V2, low_constraint$V2, na.rm=T)
t.test(is.na(high_constraint$V3), is.na(low_constraint$V3))
t.test(high_constraint$V3, low_constraint$V3, na.rm = T)

#Average breakdown of variants
x = length(high_genes$impact[high_genes$impact == "HIGH"])
x/length(high_genes$impact)
x = length(high_genes$impact[high_genes$impact == "MODERATE"])
x/length(high_genes$impact)
x = length(high_genes$impact[high_genes$impact == "LOW"])
x/length(high_genes$impact)
x = length(high_genes$impact[high_genes$impact == "MODIFIER"])
x/length(high_genes$impact)

x = length(low_genes$impact[low_genes$impact == "HIGH"])
x/length(low_genes$impact)
x = length(low_genes$impact[low_genes$impact == "MODERATE"])
x/length(low_genes$impact)
x = length(low_genes$impact[low_genes$impact == "LOW"])
x/length(low_genes$impact)
x = length(low_genes$impact[low_genes$impact == "MODIFIER"])
x/length(low_genes$impact)

#Compare impact between variant frequency regions
chisq.test(se_ann$group, se_ann$impact)

########
#######
se_e <-gb %>% group_by(consequence) %>% summarise(se_e = mean(AF,na.rm=T))
se_e <- as.data.frame(se_e)
ann_e <-gb %>% group_by(consequence_ann) %>% summarise(ann_e = mean(AF,na.rm=T))
ann_e <- as.data.frame(ann_e)
colnames(ann_e) <- c("consequence", "ann_e")
levels(ann_e$consequence) <- c("frameshift_variant", "stop_gained", "stop_lost", "splice_region_variant", "gene_fusion", "start_lost")
ann_e <- rbind(ann_e, c("splice_region_variant", "NA"))
ann_e <- rbind(ann_e, c("gene_fusion", "NA"))
ann_e <- rbind(ann_e, c("start_lost", "NA"))
se_ann <- merge(se_e, ann_e, identity = "consequence")

dfm1 <- melt(se_ann, id.vars = "consequence", na.rm = TRUE)
dfm1$value <- as.numeric(dfm1$value)

x = ggplot(dfm1, aes(x = consequence, y = value, fill = variable, color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Allele frequency") + 
  xlab("Variant consequence") + scale_x_discrete(labels=c("frameshift\nvariant", 
                                                          "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,0.25)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/gb_consequence_AF.jpeg", x, base_height = 3, base_width = 4)

x = ggplot(dfm1, aes(x = consequence, y = value, fill = variable, color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Allele frequency") + 
  xlab("Variant consequence") + scale_x_discrete(labels=c("frameshift\nvariant", 
                                                          "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,0.25)) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/gb_consequence_AF.jpeg", x, base_height = 3, base_width = 4)

#Get scatter plot 
x = ggplot(gb, aes(x = consequence, y = consequence_ann, color= group)) +
  geom_jitter(aes(color=group)) + ylab("Consequence (ANNOVAR)") + 
  xlab("Consequence (SnpEff)") + scale_x_discrete(labels=c("frameshift\nvariant", 
                                                           "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_discrete(labels=c("frameshift\nvariant", "stop\ngained", "stop\nlost")) +
  scale_color_manual(values =c("darkgrey", "blue", "hotpink"))+ 
  theme(panel.background = element_blank(), 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/gb_consequence_ANN_SE.jpeg", x, base_height = 3, base_width = 6)


un <- read.table("unique_gb_brief.txt", header=T,sep="\t")

length(un$chrom)
mean(un$AF)
range(un$AF)

table(un$Breed)
table(un$SNP)

un_hist <- ggplot(un, aes(x = Breed)) +
  geom_histogram(stat = "count") + ylab("Total unique\ngb variants") + 
  xlab("Breed") + scale_y_continuous(limits = c(0,250)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/unique_gb_breed.jpeg", un_hist, base_height = 2, base_width = 4.5)

my_pal <- c("dodgerblue2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#FF3300", "#33FFFF", 
            "#CC79A7", "#000000","hotpink","#CC0033")

un$fill <- un$Breed
levels(un$fill) <- c("dodgerblue2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#FF3300", "#33FFFF", 
                     "#CC79A7", "#000000","hotpink","#CC0033")

#un_hist <- ggplot(un, aes(x = Breed, fill = group)) +
un_hist <- ggplot(un, aes(x = Breed)) +
  geom_histogram(stat = "count", aes(fill = Breed)) + ylab("Total unique\nLOF variants") + 
  xlab("Breed") + scale_y_continuous(limits = c(0,250)) + 
  theme(panel.background = element_blank(), plot.background = element_blank(),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/unique_LOF_breed.jpeg", un_hist, base_height = 2, base_width = 4.5)



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


