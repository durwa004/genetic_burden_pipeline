#Only include autosomes and chr X (not MT and unplaced contigs)
###bcftools
bcftools <- read.table("bcftools_number_of_variants.txt", header=T)
bcf_snp <- sum(bcftools$no_SNPs)
bcf_indel <- sum(bcftools$no_indels)
bcf_tstv <- mean(bcftools$tstv)

###gatk
gatk <- read.table("gatk_number_of_variants.txt", header=T)
gatk_snp <- sum(gatk$no_SNPs)
gatk_indel <- sum(gatk$no_indels)
gatk_tstv <- mean(gatk$tstv)

###union
union <- read.table("union_number_of_variants.txt", header=T)
union_snp <- sum(union$no_SNPs)
union_indel <- sum(union$no_indels)
union_tstv <- mean(union$tstv)

###intersect
intersect <- read.table("intersect_number_of_variants.txt", header=T)
intersect_snp <- sum(intersect$no_SNPs)
intersect_indel <- sum(intersect$no_indels)
intersect_tstv <- mean(intersect$tstv)

#Create venn diagram
source("http://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph"))
library(devtools)
install_github("js229/Vennerable")
library(Vennerable)

#SNPs
vcombo <- Venn(SetNames=c("gatk","bcftools"),Weight=c(0,gatk_snp-intersect_snp,bcf_snp-intersect_snp,intersect_snp))
plot(vcombo)

##This works to remove text:

jpeg("HC_bcftools_intersect_venn.jpeg",width=6,height=6,units="in",res=1350)
plot(vcombo,show=list(SetLabels=FALSE,FaceText=FALSE,Faces=FALSE))
dev.off()

#indels
vcombo <- Venn(SetNames=c("gatk","bcftools"),Weight=c(0,gatk_indel-intersect_indel,bcf_indel-intersect_indel,intersect_indel))
plot(vcombo)

##This works to remove text:
jpeg("HC_bcftools_intersect_indel_venn.jpeg",width=6,height=6,units="in",res=1350)
plot(vcombo,show=list(SetLabels=FALSE,FaceText=FALSE,Faces=FALSE))
dev.off()