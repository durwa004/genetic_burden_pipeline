library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/")

#GB
data = read.table("ann_se_gb_by_individual.txt", header=F) # V3 = het, V4 = hom, V5 = missing?
data$total = data$V3 + data$V4
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
DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "missing", "total")
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
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#Number of homozygous variants
x1 <- plot(n_hom_gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of number of homozygous genetic burden") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

x_c <- plot_grid(x,x1,labels = "AUTO", ncol = 1)
save_plot("../Paper_2019/Nature_genetics/Figures/gb_gb_hom_EMMEANS.tiff", x_c, base_height = 7, base_width = 8)


#LOF EMMEANS
data = read.table("ann_se_lof_by_individual.txt", header=F)
data$total = data$V3 + data$V4
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


DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "missing", "total")
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
x <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#Number of homozygous variants
x1 <- plot(n_hom_gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of number of homozygous LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
x_c <- plot_grid(x,x1,labels = "AUTO", ncol = 1)
save_plot("../Paper_2019/Nature_genetics/Figures/lof_lof_hom_EMMEANS.tiff", x_c, base_height = 7, base_width = 8)



###Without accounting for DOC ###.
x = ggplot(data_br, aes(x=V2, y=total)) + theme_bw() + ylab("Genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V4)) + theme_bw() + ylab("Homozygous genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("homozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V3)) + theme_bw() + ylab("Heterozygous genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("heterozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)


  
##LOF without DOC


x = ggplot(data_br, aes(x=V2, y=total)) + theme_bw() + ylab("Genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V4)) + theme_bw() + ylab("Homozygous genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("homozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)

x = ggplot(data_br, aes(x=V2, y=V3)) + theme_bw() + ylab("Heterozygous genetic burden") + 
  xlab("Breed") + geom_boxplot() +scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("heterozygous_GB_by_breed.tiff", x, base_height = 3.5, base_width = 6)
