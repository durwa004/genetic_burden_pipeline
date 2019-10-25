library(dplyr)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/")

###Common breed and rare population variants
common_breed <- read.table("common_breed_rare_pop_snpeff_info.txt", header = T, quotes = FALSE)
                                            
#Impact
cb_t <- table(common_breed$impact)
cb_df <- as.data.frame(cb_t)
cb_df$Var1 <- factor(cb_df$Var,levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
bp <- ggplot(cb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Impact") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
  panel.border = element_blank(), axis.line.x = element_line(),
  axis.line.y = element_line(), axis.text.x = element_text(), 
  axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("common_breed_rare_pop_variant_impact.tiff", bp, base_height = 6, base_width = 6)

#Variant effect
cb_t <- table(common_breed$consequence)
cb_df <- as.data.frame(cb_t)
bp <- ggplot(cb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Variant effect") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
                                          panel.border = element_blank(), axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
                                          axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("common_breed_rare_pop_variant_effect.tiff", bp, base_height = 6, base_width = 6)

#Number of genes - do this in python
length(unique(common_breed$gene)) #1478
common_breed[which(common_breed$impact == "HIGH"),]
common_breed[which(common_breed$impact == "MODERATE"),]
common_breed[which(common_breed$impact == "LOW"),]
table(common_breed$breed_common[which(common_breed$impact == "LOW")])


###Rare breed and Common population variants
rare_breed <- read.table("rare_breed_common_pop_snpeff_info.txt", header = T)

#Impact
rb_t <- table(rare_breed$impact)
rb_df <- as.data.frame(rb_t)
rb_df$Var1 <- factor(rb_df$Var,levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
bp <- ggplot(rb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Impact") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
                                          panel.border = element_blank(), axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), 
                                          axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("rare_breed_common_pop_variant_impact.tiff", bp, base_height = 6, base_width = 6)

#Variant effect
rb_t <- table(rare_breed$consequence)
rb_df <- as.data.frame(rb_t)
bp <- ggplot(rb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Variant effect") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
                                          panel.border = element_blank(), axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
                                          axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("rare_breed_common_pop_variant_effect.tiff", bp, base_height = 6, base_width = 6)

#Number of genes
length(unique(rare_breed$gene)) #26,007

#Breed distribution
rb_b <- table(rare_breed$breed_rare)
rbb_df <- as.data.frame(rb_b)
bp <- ggplot(rbb_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Breeds with rare variants compared to population") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
                                          panel.border = element_blank(), axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
                                          axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("rare_breed_common_pop_breed.tiff", bp, base_height = 6, base_width = 6)

#Number of coding variants
length(unique(rare_breed$gene[which(rare_breed$impact != "MODIFIER")]))
length(rare_breed$gene[which(rare_breed$impact != "MODIFIER")])
rb_e <- table(rare_breed[which(rare_breed$impact != "MODIFIER")])
rbe_df <- as.data.frame(rare_breed[which(rare_breed$impact != "MODIFIER"),])
rbe_t <- table(rbe_df$consequence)
rbet_df <- as.data.frame(rbe_t)
bp <- ggplot(rbet_df,aes(x=Var1,y=Freq)) + theme_bw() + 
  ylab("Frequency") + xlab("Consequence of rare coding variants in breed and common in population") + geom_bar(stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank(),
                                          panel.border = element_blank(), axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(angle=90), 
                                          axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("rare_breed_common_pop_coding_consequence.tiff", bp, base_height = 6, base_width = 6)


table(rare_breed$breed_rare[which(rare_breed$impact == "HIGH")])
rare_breed[which(rare_breed$impact == "HIGH"),]
table(rare_breed$breed_rare[which(rare_breed$impact == "MODERATE")])
rare_breed[which(rare_breed$impact == "MODERATE"),]
rare_breed[which(rare_breed$impact == "LOW"),]
table(common_breed$breed_rare[which(common_breed$impact == "LOW")])
