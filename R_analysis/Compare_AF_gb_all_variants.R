setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/")

#Import data
gb <- read.table("AF_gb_variants.txt", header=T)
gb$type <- rep("GB",length(gb))
variants <- read.table("AF_all_variants.txt", header=T)
variants$type <- rep("all",length(variants))

gb_variants <- merge(gb,variants, by="")
#Compare the means
t.test(gb$GB_AF, variants$all_AF)
range(gb$GB_AF)
range(variants$all_AF)

