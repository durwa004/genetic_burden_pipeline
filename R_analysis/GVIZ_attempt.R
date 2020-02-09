###Create Figure 3: This code creates the HapQTL figures aligned with the equine genome
library("Gviz")
library("biomaRt")
library(karyoploteR)
library(cowplot)
setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/")
lof <- read.table("lof_details_brief.txt", header=T,sep="\t")

#helpful links for Gvis
#https://www.rdocumentation.org/packages/Gviz/versions/1.16.3/topics/plotTracks
#https://rdrr.io/bioc/Gviz/man/settings.html

#This is all part of Sam's code to lode in the genome information required to use the genome browser tools.
# This line sets up the EquCab2 chromosome sizes
#make GRanges object for EC3
equCab3 <- makeGRangesFromDataFrame(data.frame(chr=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31","chrX"),start=c(rep(0,32)),stop=c(188260000,121350000,121350000,109460000,96760000,87230000,100790000,97560000,85790000,85160000,61680000,36990000,43780000,94600000,92850000,88960000,80720000,82640000,62680000,65340000,58980000,50930000,55560000,48290000,40280000,43150000,40250000,47350000,34780000,31400000,26000000,128210000)))

gtrack <- GenomeAxisTrack(fontfamily="Arial",fontfamily.title="Arial",col=1,
                          fontsize=8,fontcolor="black",lwd=1)

#Plot whole karyoplot
#load cytoband information
ideo <- read.table("../horse_bands.txt",header=T) 
#plot karyotype; cex=0.75 is for smaller chromosome names
kp <- plotKaryotype(genome=equCab3,cytobands=ideo,chromosomes="all",use.cache=F,cex=0.75)
levels(lof$chrom) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31","chrX")
lof$Snps <- paste(lof$chrom, ":", lof$pos, "-", lof$pos,sep="")

SNPs <- GRanges(lof$Snps)
# plot inversion ranges on the ideogram in orange
x <- kpPlotRegions(kp,SNPs, col = "orange")
save_plot("../../../Abstracts/PAG_2020/Presentation/LOF_per_chromosome.jpeg", x, base_height = 8, base_width = 3)


#Get region and chromosomes
library(Gviz)
library(biomaRt)

#helpful links for Gvis
#https://www.rdocumentation.org/packages/Gviz/versions/1.16.3/topics/plotTracks
#https://rdrr.io/bioc/Gviz/man/settings.html

equCab3 <- makeGRangesFromDataFrame(data.frame(chr=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chr30","chr31","chrX"),start=c(rep(0,32)),stop=c(188260000,121350000,121350000,109460000,96760000,87230000,100790000,97560000,85790000,85160000,61680000,36990000,43780000,94600000,92850000,88960000,80720000,82640000,62680000,65340000,58980000,50930000,55560000,48290000,40280000,43150000,40250000,47350000,34780000,31400000,26000000,128210000)))

gtrack <- GenomeAxisTrack(fontfamily="Arial",fontfamily.title="Arial",col=1,
                          fontsize=8,fontcolor="black",lwd=1)

#plotTracks(gtrack,from=80500000,to=81809066)

setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/")
lof <- read.table("lof_details_brief.txt", header=T,sep="\t")

#This part of the script uses cordinates from NCBI manually put into a dataframe
Gene_Cordinates_DF <- read.table("ProteinTable145_358900.txt",header=F, sep = "\t")
Gene_Cordinates_DF$V1 <- rep("chr20",1864)
colnames(Gene_Cordinates_DF)<-c("chrom", "chrom_rs", "start", "stop", "strand", "Geneid", 
                                "locus", "protein", "length", "prot_name")
Gene_Cordinates_DF <- droplevels(Gene_Cordinates_DF %>% filter(Gene_Cordinates_DF$start > 25000000 & Gene_Cordinates_DF$stop < 34500000))

Gene_Cordinates <- makeGRangesFromDataFrame(Gene_Cordinates_DF,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                            seqnames.field="chrom",start.field="start",end.field="stop")

aTrack <- AnnotationTrack(Gene_Cordinates,name="NCBI EquCab3.0",id=Gene_Cordinates$locus)

aTrack.groups <- AnnotationTrack(Gene_Cordinates,name="NCBI EquCab3.0",id=Gene_Cordinates$locus,
                                 group=Gene_Cordinates$locus,stackHeight=0.25,fontcolor="black",fontcolor.group="black",
                                 fontfamily.group="Arial",fontsize=8,fontsize.group=10,shape="box",groupAnnotation="group",fill="black",
                                 col.title="black",col="black",background.title="transparent",fontface.title=1)
plotTracks(aTrack)

lof_mult <- droplevels(y %>% filter(count >10))
lof_m <- droplevels(lof %>% filter(chrom == "NC_009163.3") %>% filter(gene_ann %in% lof_mult$Gene))
lof_m$chr <- rep("chr20", 81)

lof_m1 <- droplevels(lof_m[,c(2,3,7)])
lof_m1$chrom <- rep("chr20", 81)
data2 <- makeGRangesFromDataFrame(lof_m1,keep.extra.columns=TRUE,ignore.strand=TRUE,
                                  seqnames.field="chrom",start.field="pos",end.field="pos")


dTrack2 <- DataTrack(data2,name="Allele frequency",cex=0.25,col.baseline=c(1,"red"),ylim = c(0,0.8),
                     lwd.baseline=c(1,2),from=25000000, to=34500000,fontcolor="black",
                     col.symbol=1,col.axis=1,cex.axis=1,fontfamily="Arial",fontfamily.title="Arial",
                     fontsize=8,col.title="black",background.title="transparent",fontface.title=1)

plotTracks(dTrack2)
plotTracks(list(dTrack2,aTrack,gtrack),fontsize=8,innerMargin =-1,title.width = 0.75,from=24800000, to=35000000)

#Get region around 25000000
Gene_Cordinates_DF <- read.table("ProteinTable145_358900.txt",header=F, sep = "\t")
Gene_Cordinates_DF$V1 <- rep("chr20",1864)
colnames(Gene_Cordinates_DF)<-c("chrom", "chrom_rs", "start", "stop", "strand", "Geneid", 
                                "locus", "protein", "length", "prot_name")
Gene_Cordinates_DF <- droplevels(Gene_Cordinates_DF %>% filter(Gene_Cordinates_DF$start > 25000000 & Gene_Cordinates_DF$stop < 26000000))

Gene_Cordinates <- makeGRangesFromDataFrame(Gene_Cordinates_DF,keep.extra.columns=TRUE,ignore.strand=FALSE,
                                            seqnames.field="chrom",start.field="start",end.field="stop")

aTrack <- AnnotationTrack(Gene_Cordinates,name="NCBI EquCab3.0",id=Gene_Cordinates$locus)

aTrack.groups <- AnnotationTrack(Gene_Cordinates,name="NCBI EquCab3.0",id=Gene_Cordinates$locus,
                                 group=Gene_Cordinates$locus,stackHeight=0.25,fontcolor="black",fontcolor.group="black",
                                 fontfamily.group="Arial",fontsize=8,fontsize.group=10,shape="box",groupAnnotation="group",fill="black",
                                 col.title="black",col="black",background.title="transparent",fontface.title=1)
plotTracks(aTrack.groups)

dTrack2 <- DataTrack(data2,name="Allele frequency",cex=0.25,col.baseline=c(1,"red"),ylim = c(0,0.8),
                     lwd.baseline=c(1,2),from=25000000, to=26000000,fontcolor="black",
                     col.symbol=1,col.axis=1,cex.axis=1,fontfamily="Arial",fontfamily.title="Arial",
                     fontsize=8,col.title="black",background.title="transparent",fontface.title=1)

plotTracks(list(dTrack2,aTrack.groups,gtrack),fontsize=8,innerMargin =-1,title.width = 0.75,from=24800000, to=26000000)


