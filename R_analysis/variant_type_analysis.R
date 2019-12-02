library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/variant_type_analysis/")

########Annovar - remove chromosome unknown variants
annovar <- read.table("annovar_variant_type_all.txt", header=F)
ann_t <- table(annovar$V5)
ann_df <- as.data.frame(ann_t)
bp <- ggplot(ann_df,aes(x="",y=Freq,fill=Var1)) + geom_bar(width=1,stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank())
ann_pie <- bp + coord_polar("y",start=0)
save_plot("annovar_variant_type.tiff", ann_pie, base_height = 6, base_width = 6)

#Variant effect by mean AF
x = ggplot(annovar, aes(x = V5, y = V4)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Variant effect") + geom_bar(stat = "summary", fun.y = "mean") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("annovar_variant_effect_by_AF.tiff", x, base_height = 3.5, base_width = 6)

#Impact by mean AF
#This doesn't work - because the exonic variants aren't annotated, therefore all variants are unknown
#annovar$V6 <- factor(annovar$V6, levels = c("HIGH", "MODERATE", "LOW", "UNKNOWN"))
#x = ggplot(ann_t, aes(x = V6, y = V4)) + theme_bw() + ylab("Mean allele frequency") + 
#xlab("Impact") + geom_bar(stat = "summary", fun.y = "mean") + 
 # theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
  #      axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
#Create my own df (using the high/mod/low AF from the coding variant output
ann_i <- c("HIGH", "MODERATE", "LOW", "UNKNOWN")
ann_i <- factor(ann_i, levels = c("HIGH", "MODERATE", "LOW", "UNKNOWN"))
mean(annovar_coding$V4[annovar_coding$V6 == "HIGH"])
mean(annovar_coding$V4[annovar_coding$V6 == "LOW"])
mean(annovar$V4)
ann_i_af <- c(0.106504, 0.0681152, 0.03845399, 0.1290549)
ann_df <- data.frame(ann_i, ann_i_af)
x = ggplot(ann_df, aes(x = ann_i, y = ann_i_af)) + theme_bw() + ylab("Mean allele frequency") + 
xlab("Impact") + geom_bar(stat = "summary", fun.y = "mean") + 
 theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
      axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("annovar_impact_by_AF_all.tiff", x, base_height = 3.5, base_width = 6)



########SnpEff - remove chromosome unknown variants
#Need to recode the variant type in python before using this(too many options)
SnpEff <- read.table("SnpEff_variant_type_all.txt", header=F,sep="\t")
SE_t <- table(SnpEff$V5)
SE_df <- as.data.frame(SE_t)
bp <- ggplot(SE_df,aes(x="",y=Freq,fill=Var1)) + geom_bar(width=1,stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank())
SE_pie <- bp + coord_polar("y",start=0)
save_plot("SnpEff_variant_type.tiff", SE_pie, base_height = 6, base_width = 6)

#Variant effect by mean AF
se_e <-SnpEff %>% group_by(V5) %>% summarise(se_e = mean(V4,na.rm=T))
se_e <- as.data.frame(se_e)
x = ggplot(se_e, aes(x = V5, y = se_e)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Variant effect") + geom_bar(stat = "identity") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("SnpEff_variant_effect_by_AF.tiff", x, base_height = 3.5, base_width = 6)

#Impact by mean AF
se_i <-SnpEff %>% group_by(V6) %>% summarise(se_i = mean(V4,na.rm=T))
se_i <-as.data.frame(se_i)
se_i$V6 <- factor(se_i$V6,levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))

x = ggplot(se_i, aes(x = V6, y = se_i)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Impact") + geom_bar(stat = "identity") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("SnpEff_impact_by_AF_all.tiff", x, base_height = 3.5, base_width = 6)


########Annovar coding
annovar_coding <- read.table("annovar_variant_type_coding.txt", header=F,sep="\t")
ann_t <- table(annovar_coding$V5)
ann_df <- as.data.frame(ann_t)
bp <- ggplot(ann_df,aes(x="",y=Freq,fill=Var1)) + geom_bar(width=1,stat = "identity") + 
  scale_y_continuous(labels=comma) + theme(legend.title = element_blank())
ann_pie <- bp + coord_polar("y",start=0)
save_plot("annovar_coding_variant_type.tiff", ann_pie, base_height = 6, base_width = 6)

#Coding variants impact by mean AF
annovar_coding$V6 <- factor(annovar_coding$V6,levels = c("HIGH", "MODERATE", "LOW", "UNKNOWN"))
x = ggplot(annovar_coding, aes(x = V6, y = V4)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Impact") + geom_bar(stat = "summary", fun.y = "mean") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("annovar_impact_by_AF.tiff", x, base_height = 3.5, base_width = 6)

#Variant effect by mean AF
x = ggplot(annovar_coding, aes(x = V5, y = V4)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Variant effect") + geom_bar(stat = "summary", fun.y = "mean") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("annovar_coding_variant_effect_by_AF.tiff", x, base_height = 3.5, base_width = 6)

########Snpeff coding
SnpEff_coding <- read.table("SnpEff_variant_type_coding.txt", header=F,sep="\t")
SE_t <- table(SnpEff_coding$V5)
SE_df <- as.data.frame(SE_t)
bp <- ggplot(SE_df,aes(x="",y=Freq,fill=Var1)) + geom_bar(width=1,stat = "identity") + 
  scale_y_continuous(labels=comma)+ theme(legend.title = element_blank())
SE_pie <- bp + coord_polar("y",start=0)
save_plot("SnpEff_variant_type_coding.tiff", SE_pie, base_height = 6, base_width = 6)

#Coding variants impact by mean AF
se_i <-SnpEff_coding %>% group_by(V6) %>% summarise(se_i = mean(V4,na.rm=T))
se_i <- as.data.frame(se_i)
se_i$V6 <- factor(se_i$V6,levels = c("HIGH", "MODERATE", "LOW"))
x = ggplot(se_i, aes(x = V6, y = se_i)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Impact") + geom_bar(stat = "identity") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("SnpEff_impact_by_AF_coding.tiff", x, base_height = 3.5, base_width = 6)

#Variant effect by mean AF
se_e <-SnpEff_coding %>% group_by(V5) %>% summarise(se_e = mean(V4,na.rm=T))
se_e <- as.data.frame(se_e)
x = ggplot(se_e, aes(x = V5, y = se_e)) + theme_bw() + ylab("Mean allele frequency") + 
  xlab("Variant effect") + geom_bar(stat = "identity") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("SnpEff_coding_variant_effect_by_AF.tiff", x, base_height = 3.5, base_width = 6)

