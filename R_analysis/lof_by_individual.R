library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper/")

#GB
data = read.table("lof_by_individual.txt", header=F) # V3 = het, V4 = hom, V5 = missing?
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
x3 <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits = c(540,900)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"), plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Number of homozygous variants
x4 <- plot(n_hom_gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of homozygous LOF variants") + 
  ylab("Breed") +scale_x_continuous(labels=comma, limits =c(100,220)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"), plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Combine with GB_by_individual.R for nature paper
x_c <- plot_grid(x,x1,labels = "AUTO", ncol = 1)
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_EMMEANS.jpeg", x_c, base_height = 7, base_width = 8)

#Plot bar plot of LOF variants
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}

#Need to categorize DOC
gb_br$DOC <- cut(gb_br$nuclear_placed_DOC, 
                       breaks = c(-Inf, 4, 8, 12, 16, 20, 30, 40, Inf), 
                       labels = c("<4", "4-7", "8-11", "12-15", "16-19","20-29","30-39",">40"), 
                       right = FALSE)
gb_hist <- ggplot(gb_br, aes(x = total,fill = DOC,color = DOC)) +
  geom_histogram(position = "identity", binwidth = 10) + facet_grid(scales = "free_y",breed ~ .) + 
  ylab("Frequency") + xlab("Number of LOF variants") + scale_y_continuous(labels = formatter(nsmall = 1)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
      plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
      legend.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
      legend.position = "bottom", panel.grid = element_blank(), panel.border = element_blank(), 
      axis.line.x = element_line(), axis.line.y = element_line(), 
      axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
      axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/LOF_histogram.jpeg", gb_hist, base_height = 10, base_width = 4)


#GET INFORMATION RE Ne

#Get values from emmeans
gb_emm
n_hom_gb_emm

gb_t <- as.data.frame(table(gb_br$breed))
colnames(gb_t) <- c("breed", "nvariants")
gb_t$nvariants <- c(755,807,752,817,615,"NA",744,706,709,626,729)
gb_t$nvariants <- as.numeric(gb_t$nvariants)
gb_t$nhom <- c(180,192,199,199,127,"NA",159,161,155,128,165)
gb_t$nhom <- as.numeric(gb_t$nhom)

# Add in Ne from Jessica's paper
gb_t$ne1 <- "NA"
gb_t$ne1[gb_t$breed == "Arabian"] <- 346
gb_t$ne1[gb_t$breed == "Belgian"] <- 431
gb_t$ne1[gb_t$breed == "Icelandic"] <- 555
gb_t$ne1[gb_t$breed == "Morgan"] <- 448
gb_t$ne1[gb_t$breed == "QH"] <- 426
gb_t$ne1[gb_t$breed == "STB"] <- 290
gb_t$ne1[gb_t$breed == "TB"] <- 190
gb_t$ne1 <- as.numeric(gb_t$ne1)

cor.test(x=gb_t$ne1,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne1,y=gb_t$nhom, method = "pearson", use='complete.obs')

#Add in effective population size (from Sam's paper)
gb_t$ne <- "NA"
gb_t$ne[gb_t$breed == "Arabian"] <- 3561
gb_t$ne[gb_t$breed == "Belgian"] <- 3570
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "Morgan"] <- 4481
gb_t$ne[gb_t$breed == "Icelandic"] <- 2736
gb_t$ne[gb_t$breed == "QH"] <- 6516
gb_t$ne[gb_t$breed == "STB"] <- 2528
gb_t$ne[gb_t$breed == "TB"] <- 1784
gb_t$ne[gb_t$breed == "WP"] <- 5625

gb_t$ne <- as.numeric(gb_t$ne)
cor.test(x=gb_t$ne,y=gb_t$nvariants, method = "pearson", use='complete.obs')
cor.test(x=gb_t$ne,y=gb_t$nhom, method = "pearson", use='complete.obs')


####Plot EMMEANS
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
#Number of variants on same axis (poster)
labels <- c("1" = "Hom. LOF variants", "2" = "All LOF variants")
x1 <- plot(gb_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of LOF variants") + 
  scale_fill_manual(values = c(rgb(155/255,218/255,233/255,1), rgb(122/255,0/255,25/255,1))) + 
  ylab("Breed") +scale_x_continuous(labels=comma) + 
  facet_grid(v_type~., labeller=labeller(v_type = labels)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
save_plot("../../../Abstracts/PAG_2020/Poster/LOF_nhom_EMMEANS.tiff", x1, base_height = 3.5, base_width = 5)



#Number of variants - with different axis sizes (presentation)
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
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_EMMEANS.tiff", x_c, base_height = 6, base_width = 6)



#Get type of variants
#Have to go through and delete the # in excel
lof <- read.table("lof_details_brief.txt", header=T,sep="\t")

table(lof$group)

mean(lof$AF)
range(lof$AF)

table(lof$SNP)
table(lof$SNP_ann)

table(lof$consequence)
table(lof$consequence_ann)

y <- as.data.frame(table(lof$consequence))
colnames(y) <- c("consequence", "SnpEff")
z <- as.data.frame(table(lof$consequence_ann))
colnames(z) <- c("consequence", "ANNOVAR")
levels(z$consequence) <- c("frameshift_variant", "stop_gained", "stop_lost", "splice_region_variant", "gene_fusion", "start_lost")
z <- rbind(z, c("splice_region_variant", "NA"))
z <- rbind(z, c("gene_fusion", "NA"))
z <- rbind(z, c("start_lost", "NA"))

se_ann1 <- merge(y, z, identity = "consequence")
dfm <- melt(se_ann1,  id.vars = "consequence", na.rm = TRUE)
dfm$value <- as.numeric(dfm$value)

lof_bar <- ggplot(dfm, aes(x = consequence, y = value, fill = variable,color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Frequency") + 
  xlab("Number of LOF variants") + scale_x_discrete(labels=c("frameshift\nvariant", 
  "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,3000)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/LOF_consequence.jpeg", lof_bar, base_height = 3, base_width = 4)

lof_bar <- ggplot(dfm, aes(x = consequence, y = value, fill = variable,color = variable)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.8, alpha = 1) + ylab("Frequency") + 
  xlab("Number of LOF variants") + scale_x_discrete(labels=c("frameshift\nvariant", 
  "gene\nfusion", "splice region\nvariant", "start\nlost", "stop\ngained", "stop\nlost")) +
  scale_y_continuous(limits = c(0,3000)) +
  theme(panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_consequence.jpeg", lof_bar, base_height = 3, base_width = 6)



se_e <-lof %>% group_by(consequence) %>% summarise(se_e = mean(AF,na.rm=T))
se_e <- as.data.frame(se_e)
ann_e <-lof %>% group_by(consequence_ann) %>% summarise(ann_e = mean(AF,na.rm=T))
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
save_plot("../../../Abstracts/PAG_2020/Poster/LOF_consequence_AF.jpeg", x, base_height = 3, base_width = 4)

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
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_consequence_AF.jpeg", x, base_height = 3, base_width = 4)

#Get scatter plot 
x = ggplot(lof, aes(x = consequence, y = consequence_ann, color= group)) +
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
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_consequence_ANN_SE.jpeg", x, base_height = 3, base_width = 6)


un <- read.table("unique_lof_brief.txt", header=T,sep="\t")

length(un$chrom)
mean(un$AF)
range(un$AF)

table(un$Breed)
table(un$SNP)

un_hist <- ggplot(un, aes(x = Breed)) +
  geom_histogram(stat = "count") + ylab("Total unique\nLOF variants") + 
  xlab("Breed") + scale_y_continuous(limits = c(0,250)) +
  theme(panel.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"), 
        plot.background = element_rect(fill=rgb(253/255,244/255,210/255,1), colour = "NA"),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.ticks.x = element_blank(), axis.title = element_text(size=12,face="bold"),
        axis.text.x = element_text(angle = 90), plot.margin=unit(c(0,0,0,0), "null"))
save_plot("../../../Abstracts/PAG_2020/Poster/unique_LOF_breed.jpeg", un_hist, base_height = 2, base_width = 4.5)

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



