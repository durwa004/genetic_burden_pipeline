library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)
library(devtools)
library(ggpubr)
library(emmeans)
library(scales)
library(dvmisc)
library(tidyr)
library(gt)
library(tidyverse)

setwd("/Users/durwa004/Documents/Research/GB_project")
####Figure out differences in the number of variants per breed
###Individual GB
#use ind


#Look for association between breed and genetic burden accounting for DOC
#Add in DOC and breed info 
bcftools <- read.table("ind_number_of_variants.txt", header=T)

#Add in breed information
breed <- read.table("horse_genomes_breeds_all.txt", header = T)
colnames(breed) <- c("ID", "breed")
bcftools <- merge(bcftools, breed, "ID")
bcftools <- bcftools %>%
  mutate(across('Breed', str_replace_all, "_", " ")) %>%
  mutate(across('Breed', str_replace, "_", " ")) %>%
  mutate(across('Breed', str_replace, "_", " ")) %>%
  mutate(across('Breed', str_replace, "_x", " x")) %>%
  mutate(across('Breed', str_replace, "Uknown", "Unknown")) %>%
  mutate(across('Breed', str_replace, "Trakenher", "Trakehner")) %>%
  mutate(across('Breed', str_replace, "Halflinger", "Haflinger")) %>%
  mutate(across('Breed', str_replace, "Hanovarian", "Hanoverian")) %>%
  mutate(across('Breed', str_replace, "ColdBlood", "Coldblood")) %>%
  mutate(across('Breed', str_replace, "Other horse", "Unknown")) %>%
  mutate(across('Breed', str_replace, "Pony", "Unknown")) %>%
  mutate(across('Breed', str_replace, "QH", "Quarter Horse")) %>%
  mutate(across('Breed', str_replace, "STB", "Standardbred")) %>%
  mutate(across('Breed', str_replace, "TB", "Thoroughbred")) %>%
  mutate(across('Breed', str_replace, "WarmBlood", "Warmblood")) %>%
  mutate(across('Breed', str_replace, "TWH", "Tennessee Walking Horse")) %>%
  mutate(across('Breed', str_replace, "duelmener", "Dulmener"))


colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "missing", "total")
gb_doc <- merge(data,DOC, by="Sample")
gb_br <- gb_doc %>% 
  filter(!grepl('Other', breed))

fit1 <- (lm(total ~ breed, data=gb_br))
fit2 <- (lm(total ~ breed + nuclear_placed_DOC, data=gb_br))
anova(fit1,fit2)

gb_hom <- gb_br[c(1,2,4,8)]
colnames(gb_hom) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")
gb_het <- gb_br[c(1,2,3,8)]
colnames(gb_het) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")
gb_nv <- gb_br[c(1,2,6,8)]
colnames(gb_nv) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")

gb_stats <- bind_rows(gb_hom, gb_nv, .id = "v_type")
gb_stats$v_type <- as.factor(gb_stats$v_type)

gb_m <- (lm(variants ~ v_type + breed + nuclear_placed_DOC,data=gb_stats))

#Get EMMEANs
gb_emm <- emmeans(gb_m, specs = "breed", weights = "proportional", 
                         type = "response", by = "v_type")

####Plot EMMEANS
#Number of variants
labels <- c("1" = "Homozygous", "2" = "Genetic burden")
x1 <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of number of variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  facet_grid(v_type~., labeller=labeller(v_type = labels)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../Paper_2019/Figures/Ind_GB_nvariants_nhom_EMMEANS.tiff", x1, base_height = 3.5, base_width = 8)

#LOF EMMEANS
data = read.table("ann_se_lof_by_individual.txt", header=F)
data$total = data$V3 + data$V4
mean(data$total) 
range(data$total)
mean(data$V3)
range(data$V3)
mean(data$V4)
range(data$V4)

DOC <- read.table("../DOC/DOC_by_horse.txt", header=T)
colnames(DOC) <- c("Sample", "total_DOC", "nuclear_placed_DOC")
colnames(data) = c("Sample", "breed", "het", "hom", "missing", "total")
lof_doc <- merge(data,DOC, by="Sample")
lof_br <- lof_doc %>% 
  filter(!grepl('Other', breed))

fit1 <- (lm(total ~ breed, data=lof_br))
fit2 <- (lm(total ~ breed + nuclear_placed_DOC, data=lof_br))
anova(fit1,fit2)

lof_hom <- lof_br[c(1,2,4,8)]
colnames(lof_hom) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")
lof_het <- lof_br[c(1,2,3,8)]
colnames(lof_het) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")
lof_nv <- lof_br[c(1,2,6,8)]
colnames(lof_nv) <- c("Sample", "breed", "variants", "nuclear_placed_DOC")

lof_stats <- bind_rows(lof_hom, lof_nv, .id = "v_type")
lof_stats$v_type <- as.factor(lof_stats$v_type)

lof_m <- (lm(variants ~ v_type + breed + nuclear_placed_DOC,data=lof_stats))

#Get EMMEANs
lof_emm <- emmeans(lof_m, specs = "breed", weights = "proportional", 
                  type = "response", by = "v_type")

####Plot EMMEANS
#Number of variants
labels <- c("1" = "Homozygous LOF", "2" = "LOF genetic burden")
x1 <- plot(lof_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of number of variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  facet_grid(v_type~., labeller=labeller(v_type = labels)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../Paper_2019/Figures/Ind_lof_nvariants_nhom_EMMEANS.tiff", x1, base_height = 3.5, base_width = 8)







