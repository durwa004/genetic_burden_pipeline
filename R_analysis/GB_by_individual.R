library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper")

#GB
data = read.table("genetic_burden_by_individual.txt", header=F) # V3 = het, V4 = hom
data$total <- data$V3 + data$V4
mean(data$total) 
range(data$total)
mean(data$V3)
range(data$V3)
mean(data$V4)
range(data$V4)

mean(data$total[data$V2 == "Arabian"]) #3874

mean(data$total[data$V2 == "Belgian"]) #4053
mean(data$total[data$V2 == "Clydesdale"]) #3763
mean(data$total[data$V2 == "Icelandic"]) #4144
mean(data$total[data$V2 == "Morgan"]) #2885
mean(data$total[data$V2 == "QH"]) #3810
mean(data$total[data$V2 == "Shetland"]) #3407
mean(data$total[data$V2 == "STB"]) #3514
mean(data$total[data$V2 == "TB"]) #3033
mean(data$total[data$V2 == "WP"]) #3603

#Look for association between breed and genetic burden accounting for DOC
#Add in DOC info
library(emmeans)
DOC <- read.table("../../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "total")
gb_doc <- merge(data,DOC, by="Sample")
gb_br <- gb_doc %>% 
  filter(!grepl('Other', breed))

fit1 <- (lm(total ~ breed, data=gb_br))
fit2 <- (lm(total ~ breed + nuclear_placed_DOC, data=gb_br))
anova(fit1,fit2)

gb_m <- (lm(total ~ breed + nuclear_placed_DOC,data=gb_br))
n_hom_gb_m <- (lm(hom ~ breed + nuclear_placed_DOC,data=gb_br))
summary(gb_m)
summary(n_hom_gb_m)

#Get EMMEANs
gb_emm <- emmeans(gb_m, "breed", weights = "proportional", type = "response")
gb_emm
n_hom_gb_emm <- emmeans(n_hom_gb_m, "breed", weights = "proportional", type = "response")
n_hom_gb_emm

####Plot EMMEANS
#Number of variants
x <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of genetic burden") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits=c(700,1050)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#Number of homozygous variants
x1 <- plot(n_hom_gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of homozygous genetic burden") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))


first_col <- plot_grid(x,x1, labels = c("A", "B"),rel_widths = 1,rel_heights = 1, ncol = 1)
second_col <- plot_grid(x3,x4, labels = c("C", "D"), rel_widths = 1,rel_heights = 1, ncol = 1)
x_c <- plot_grid(first_col,second_col, ncol = 1, rel_widths = 1, rel_heights = 1)

save_plot("../../Paper_2019/Nature_genetics/Figures/gb_gb_hom_EMMEANS.tiff", x_c, base_height = 12, base_width = 6)

#########################################

#Get type of variants
#Have to go through and delete the # using sed
gb <- read.table("genetic_burden_details_brief.txt", header=T,sep="\t")

table(gb$group)

mean(gb$AF)
range(gb$AF)

table(gb$SNP)
table(gb$SNP_ann)

table(gb$consequence)
table(gb$consequence_ann)

y <- as.data.frame(table(gb$consequence))
colnames(y) <- c("consequence", "SnpEff")
z <- as.data.frame(table(gb$consequence_ann))
colnames(z) <- c("consequence", "ANNOVAR")
levels(z$consequence) <- c("frameshift_variant", "stop_gained", "stop_lost", "splice_region_variant", "gene_fusion", "start_lost")
z <- rbind(z, c("splice_region_variant", "NA"))
z <- rbind(z, c("gene_fusion", "NA"))
z <- rbind(z, c("start_lost", "NA"))

se_ann1 <- merge(y, z, identity = "consequence")
dfm <- melt(se_ann1,  id.vars = "consequence", na.rm = TRUE)
dfm$value <- as.numeric(dfm$value)

gb_bar <- ggplot(dfm, aes(x = consequence, y = value, fill = variable,color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Frequency") + 
  xlab("Number of gb variants") + scale_x_discrete(labels=c("frameshift\nvariant", 
                                                             "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,3000)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/gb_consequence.jpeg", gb_bar, base_height = 3, base_width = 4)

gb_bar <- ggplot(dfm, aes(x = consequence, y = value, fill = variable,color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Frequency") + 
  xlab("Number of gb variants") + scale_x_discrete(labels=c("frameshift\nvariant", 
                                                             "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,3000)) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/gb_consequence.jpeg", gb_bar, base_height = 3, base_width = 6)



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



