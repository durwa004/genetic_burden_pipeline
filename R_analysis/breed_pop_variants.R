library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/")

rare_breed <- read.table("rare_breed_common_pop.txt", header=T)

rare_Arab <- rare_breed %>% 
  filter(grepl('Arabian', Breed))
length(rare_Arab$Breed)
mean(rare_Arab$AF)
mean(rare_Arab$AC)
write.table(rare_Arab, "rare_Arab_common_pop.txt",quote=FALSE)

rare_Belg <- rare_breed %>% 
  filter(grepl('Belgian', Breed))
length(rare_Belg$Breed)
mean(rare_Belg$AF)
mean(rare_Belg$AC)
write.table(rare_Belg, "rare_Belgian_common_pop.txt",quote=FALSE)

rare_Clydesdale <- rare_breed %>% 
  filter(grepl('Clydesdale', Breed))
length(rare_Clydesdale$Breed)
mean(rare_Clydesdale$AF)
mean(rare_Clydesdale$AC)
write.table(rare_Clydesdale, "rare_Clydesdale_common_pop.txt",quote=FALSE)

rare_Icelandic <- rare_breed %>% 
  filter(grepl('Icelandic', Breed))
length(rare_Icelandic$Breed)
mean(rare_Icelandic$AF)
mean(rare_Icelandic$AC)
write.table(rare_Icelandic, "rare_Icelandic_common_pop.txt",quote=FALSE)

rare_Morgan <- rare_breed %>% 
  filter(grepl('Morgan', Breed))
length(rare_Morgan$Breed)
mean(rare_Morgan$AF)
mean(rare_Morgan$AC)
write.table(rare_Morgan, "rare_Morgan_common_pop.txt",quote=FALSE)

rare_QH <- rare_breed %>% 
  filter(grepl('QH', Breed))
length(rare_QH$Breed)
mean(rare_QH$AF)
mean(rare_QH$AC)
write.table(rare_QH, "rare_QH_common_pop.txt",quote=FALSE)

rare_Shetland <- rare_breed %>% 
  filter(grepl('Shetland', Breed))
length(rare_Shetland$Breed)
mean(rare_Shetland$AF)
mean(rare_Shetland$AC)
write.table(rare_Shetland, "rare_Shetland_common_pop.txt",quote=FALSE)

rare_STB <- rare_breed %>% 
  filter(grepl('STB', Breed))
length(rare_STB$Breed)
mean(rare_STB$AF)
mean(rare_STB$AC)
write.table(rare_STB, "rare_STB_common_pop.txt",quote=FALSE)

rare_TB <- rare_breed %>% 
  filter(grepl('TB', Breed))
length(rare_TB$Breed)
mean(rare_TB$AF)
mean(rare_TB$AC)
write.table(rare_TB, "rare_TB_common_pop.txt",quote=FALSE)

rare_WP <- rare_breed %>% 
  filter(grepl('WP', Breed))
length(rare_WP$Breed)
mean(rare_WP$AF)
mean(rare_WP$AC)
write.table(rare_WP, "rare_WP_common_pop.txt",quote=FALSE)

common_breed <- read.table("common_breed_rare_pop.txt", header=T)

common_Arab <- common_breed %>% 
  filter(grepl('Arabian', Breed))
common_Arab <- common_Arab %>% filter(AC>=4)
length(common_Arab$Breed)
mean(common_Arab$AF)
mean(common_Arab$AC)
write.table(common_Arab, "common_Arab_rare_pop.txt",quote=FALSE)

common_Belg <- common_breed %>% 
  filter(grepl('Belgian', Breed))
common_Belg <- common_Belg %>% filter(AC>=2)
length(common_Belg$Breed)
mean(common_Belg$AF)
mean(common_Belg$AC)
write.table(common_Belg, "common_Belg_rare_pop.txt",quote=FALSE)

common_Clydesdale <- common_breed %>% 
  filter(grepl('Clydesdale', Breed))
common_Clydesdale <- common_Clydesdale %>% filter(AC>=2)
length(common_Clydesdale$Breed)
mean(common_Clydesdale$AF)
mean(common_Clydesdale$AC)
write.table(common_Clydesdale, "common_Clydesdale_rare_pop.txt",quote=FALSE)

common_Icelandic <- common_breed %>% 
  filter(grepl('Icelandic', Breed))
common_Icelandic <- common_Icelandic %>% filter(AC>=2)
length(common_Icelandic$Breed)
mean(common_Icelandic$AF)
mean(common_Icelandic$AC)
write.table(common_Icelandic, "common_Icelandic_rare_pop.txt",quote=FALSE)

common_Morgan <- common_breed %>% 
  filter(grepl('Morgan', Breed))
common_Morgan <- common_Morgan %>% filter(AC>=2)
length(common_Morgan$Breed)
mean(common_Morgan$AF)
mean(common_Morgan$AC)
write.table(common_Morgan, "common_Morgan_rare_pop.txt",quote=FALSE)

common_QH <- common_breed %>% 
  filter(grepl('QH', Breed))
common_QH <- common_QH %>% filter(AC>=7)
length(common_QH$Breed)
mean(common_QH$AF)
mean(common_QH$AC)

common_Shetland <- common_breed %>% 
  filter(grepl('Shetland', Breed))
common_Shetland <- common_Shetland %>% filter(AC>=6)
length(common_Shetland$Breed)
mean(common_Shetland$AF)
mean(common_Shetland$AC)

common_STB <- common_breed %>% 
  filter(grepl('STB', Breed))
common_STB <- common_STB %>% filter(AC>=5)
length(common_STB$Breed)
mean(common_STB$AF)
mean(common_STB$AC)

common_TB <- common_breed %>% 
  filter(grepl('TB', Breed))
common_TB <- common_TB %>% filter(AC>=5)
length(common_TB$Breed)
mean(common_TB$AF)
mean(common_TB$AC)

common_WP <- common_breed %>% 
  filter(grepl('WP', Breed))
common_WP <- common_WP %>% filter(AC>=2)
length(common_WP$Breed)
mean(common_WP$AF)
mean(common_WP$AC)
write.table(common_WP, "common_WP_rare_pop.txt",quote=FALSE)
