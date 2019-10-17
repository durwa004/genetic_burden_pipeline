library(ggplot2)
library(cowplot)
library(scales)

annovar_high <- 9935
annovar_moderate <- 254870
annovar_low <- 501262
annovar_modifier <- (31136852 - (annovar_high + annovar_low + annovar_moderate))

SnpEff_high <- 9891
SnpEff_moderate <- 264831
SnpEff_low <- 558532
SnpEff_modifier <- (31136852 - (SnpEff_high + SnpEff_low + SnpEff_moderate))

#Can't plot modifier - too difficult to see
program <- c(rep("SnpEFF",3),rep("Annovar",3))
impact <- rep(c("High", "Moderate", "Low"),2)
value <- c(SnpEff_high,SnpEff_moderate,SnpEff_low, annovar_high, annovar_moderate, annovar_low)
data <- data.frame(program,impact,value)

x <- ggplot(data,aes(fill=factor(impact,levels=c("High","Moderate","Low")),y = value, x=program)) + 
              geom_bar(position="stack",stat="identity") + 
  theme_bw() + ylab("Number of variants") + xlab("Program") + 
  scale_y_continuous(labels=comma) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        legend.title = element_blank(),
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
save_plot("annovar_snpeff_variants.tiff", x, base_height = 3.5, base_width = 6)

