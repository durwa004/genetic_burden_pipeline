library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/")

data = read.table("ann_se_gb_by_individual.txt", header=F)
data$total = data$V3 + data$V4
mean(data$total) # 3474
range(data$total) #777 - 5262
mean(data$V3) # 2791
range(data$V3) #615 - 4725
mean(data$V4) # 683
range(data$V4) # 162 - 1125

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

data_br <- data %>% 
  filter(!grepl('Other', V2))
kruskal.test(data_br$total, data_br$V2)

x = ggplot(data_br, aes(x=V2, y=total)) + theme_bw() + ylab("Genetic burden") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V4)) + theme_bw() + ylab("Homozygous genetic burden") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("homozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V3)) + theme_bw() + ylab("Heterozygous genetic burden") + 
  xlab("Breed") + geom_boxplot(fill=rgb(122/255,0/255,25/255,1)) +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("heterozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

#Look for association between breed and genetic burden accounting for DOC
#Add in DOC info
library(emmeans)
DOC <- read.table("../bcftools_stats_output/DOC_by_horse.txt", header=T)
colnames(data) = c("Sample", "breed", "het", "hom", "total")
gb_doc <- merge(data,DOC, by="Sample")
gb_br <- gb_doc %>% 
  filter(!grepl('Other', breed))


gb_lm <- lm(total ~ breed + nuclear_placed_DOC + breed*nuclear_placed_DOC, data = gb_br)
summary(gb_lm)
gb_em <- test(emmeans(gb_lm, ~breed, weights = "proportional", type="response"))
gb_conf <- confint(emmeans(gb_lm, ~ breed, weights = "proportional", type="response"))
gb_em_conf <- merge(gb_em,gb_conf, by="breed")


x <- ggplot(gb_em_conf) + (aes(x = breed, y = emmean.x, ymin=lower.CL, ymax=upper.CL))  + 
  geom_pointrange()  + theme(legend.title = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(), axis.line.x = element_line(),
  axis.line.y = element_line(), axis.text.x = element_text(), 
  axis.text = element_text(size=18), axis.title = element_text(size=18,face="bold")) +
  coord_flip() + theme(legend.position="none") + ylab("EMMEAN of genetic burden") + 
  xlab("Breed")
save_plot("EMMEANs_breed_GB.tiff", x, base_height = 6, base_width = 12)

bp <- ggplot(cb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Impact") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ 
#source("../genetic_burden_pipeline/R_analysis/sian-code.R")
#emms1 <- modelemms2(gb_lm,gb_doc)
#emms1x <- getsig(emms1, gb_lm)
#plotemms(emms1, digits=2) + ylab("Estimated Marginal Mean of log10 TEQ residuals")


