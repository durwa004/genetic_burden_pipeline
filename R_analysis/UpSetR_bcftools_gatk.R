#https://github.com/hms-dbmi/UpSetR
#Probably more useful for intersections of >3/4 datasets = not so much for this
library('UpSetR')
#Create a dataset with the number of variants that overlap between gatk and bcftools
#Use the results from gatk_bcftools_concordance.txt
#Number of variants
variants <- seq(1,46413496)
HC <- rep(NA,46413496)
ST <- rep(NA,46413496)
HC[0:39535390] = 1
ST[0:29882273] = 1
ST[39535391:43048393] = 1
data <- as.data.frame(variants)
data$HC <- HC
data$ST <- ST
