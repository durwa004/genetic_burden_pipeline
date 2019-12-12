setwd("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/")
data = read.table("breed_rare_common_number_of_variants.txt", header =TRUE)

mean(data$no_records)
mean(data$no_SNPs)
mean(data$no_indels)
