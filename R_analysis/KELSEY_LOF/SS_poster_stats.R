setwd("/Users/durwa004/Documents/Postdoc/Projects/LOF_variants/")
#Open table of mutations
data <- read.table("Kelsey_LOF_brief.txt", header=T)

#Are frameshift mutations more likely to be deleterious than stop codon mutations (with PROVEAN?)
fs_del <- lm(PROVEAN_score ~ consequence, data = data)
summary(fs_del)
###No significant difference between PROVEAN score between fs and sg (p =0.33)

#Are frameshift mutations more likely to be false positives than stop codon mutations?
tbl <- table(data$Real_based_on_blast., data$consequence)
chisq.test(tbl)
###Statistically significant difference between false postive frequency 
#between fs and sg (p < 0.001) - fs more likely to be false + than sg

#Are rare variants more likely to be true than common (accounting for variant type)?
rare_common <- glm(Real_based_on_blast.~ VEP + consequence, data = data, family = "binomial")
summary(rare_common)
###Statistically significant difference between false positive rate between rare and 
#common variants - rare more likely to be true (p = 0.001)

#Are there more frameshifts than stop codons in the common cf rare mutations
tbl <- table(data$consequence, data$VEP)
chisq.test(tbl)
###There is a significant difference between the type of variant between common and 
#rare variants (p < 0.001) - more frameshifts in common than rare