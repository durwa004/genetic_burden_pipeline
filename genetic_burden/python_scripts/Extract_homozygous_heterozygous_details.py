from __future__ import division  
import gzip
from statistics import mean  

#This is a very long script and likely can be tidied up at some point!!
path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis"

#Need to get header from vcf file
horse_breed = {}
with gzip.open(path + "/../bcftools_stats_output/NC_001640_1.genotyped.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[0] == "#CHROM":
            for i in range(len(line)):
                horse_breed[line[i]] = "A"

#Add in breed
with open(path + "/../bcftools_stats_output/horse_genomes_breeds_tidy.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for key in horse_breed.keys():
            if line[0] == key:
                horse_breed[key] = line[1]
#0-8 are not horses
 
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol          
##############################################################################
########################    No homozygotes    ########################

#Get list of affected genes
coding = 0
with open(path + "/no_homozygotes_details.txt", "r") as input_file, open(path + "/no_homozygotes_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[7] != "MODIFIER":
            print(line[8], file = output_file)
            coding +=1

#Create accession/gene symbol dictionary
genes = {}
with open(path + "/no_homozygotes_genes_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genes.keys():
               genes[line1[0]][line1[0]] = {}
           else:
               genes[line1[0]] = {}
               genes[line1[0]][line1[0]] = {}
       else:
           if line1[1] in genes.keys():
               genes[line1[1]][line1[0]] = {}
           else:
               genes[line1[1]] = {}
               genes[line1[1]][line1[0]] = {}
print(len(genes))

count = 0
for key in genes.keys():
    for key1 in genes[key].keys():
        count +=1
print(count)

#Get consequence and AF etc of these variants (use R)
with open(path + "/no_homozygotes_details.txt", "r") as input_file, open(path + "/no_homozygotes_summary.txt", "w") as output_file:
    input_file.readline()
    print("AC\tAF\tconsequence\timpact\tgene\tlof", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene = line[8]
        if gene == "":
            gene = "NA"
        print(line[4], line[5], line[6], line[7], gene, line[11], sep = "\t", file = output_file)

##############################################################################
########################    All homozygotes    ########################

#Get list of affected genes
coding = 0
with open(path + "/all_homozygotes_details.txt", "r") as input_file, open(path + "/all_homozygotes_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[7] != "MODIFIER":
            print(line[8], file = output_file)
            coding +=1

#Create accession/gene symbol dictionary
genes = {}
with open(path + "/all_homozygotes_genes_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genes.keys():
               genes[line1[0]][line1[0]] = {}
           else:
               genes[line1[0]] = {}
               genes[line1[0]][line1[0]] = {}
       else:
           if line1[1] in genes.keys():
               genes[line1[1]][line1[0]] = {}
           else:
               genes[line1[1]] = {}
               genes[line1[1]][line1[0]] = {}
print(len(genes))

count = 0
for key in genes.keys():
    for key1 in genes[key].keys():
        count +=1
print(count)

#Get consequence and AF etc of these variants (use R)
with open(path + "/all_homozygotes_details.txt", "r") as input_file, open(path + "/all_homozygotes_summary.txt", "w") as output_file:
    input_file.readline()
    print("AC\tAF\tconsequence\timpact\tgene\tlof", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene = line[8]
        if gene == "":
            gene = "NA"
        print(line[4], line[5], line[6], line[7], gene, line[11], sep = "\t", file = output_file)


##############################################################################
########################    Present in all    ########################
#Also need to pull out the ones that are not homozygous in all individuals
#Get list of affected genes
coding = 0
with open(path + "/variants_present_in_all_details.txt", "r") as input_file, open(path + "/variants_present_in_all_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[7] != "MODIFIER":
            print(line[8], file = output_file)
            coding +=1
print(coding)
#Create accession/gene symbol dictionary
genes = {}
with open(path + "/variants_present_in_all_genes_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genes.keys():
               genes[line1[0]][line1[0]] = {}
           else:
               genes[line1[0]] = {}
               genes[line1[0]][line1[0]] = {}
       else:
           if line1[1] in genes.keys():
               genes[line1[1]][line1[0]] = {}
           else:
               genes[line1[1]] = {}
               genes[line1[1]][line1[0]] = {}
print(len(genes))

count = 0
for key in genes.keys():
    for key1 in genes[key].keys():
        count +=1
print(count)

#Get consequence and AF etc of these variants (use R)
with open(path + "/variants_present_in_all_details.txt", "r") as input_file, open(path + "/variants_present_in_all_summary.txt", "w") as output_file:
    input_file.readline()
    print("AC\tAF\tconsequence\timpact\tgene\tlof", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene = line[8]
        if gene == "":
            gene = "NA"
        print(line[4], line[5], line[6], line[7], gene, line[11], sep = "\t", file = output_file)

#Get consequence and AF etc of these variants (use R)
with open(path + "/variants_present_in_all_not_homozygous_in_all.txt", "r") as input_file, open(path + "/variants_present_in_all_not_homozygous_in_all_summary.txt", "w") as output_file:
    input_file.readline()
    print("AC\tAF\tconsequence\timpact\tgene\tlof", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene = line[8]
        if gene == "":
            gene = "NA"
        print(line[4], line[5], line[6], line[7], gene, line[11], sep = "\t", file = output_file)



#Get number of variants with AFs >=0.50
over = []
under = []
count_over = 0
count_under = 0
with open(path + "/no_homozygotes_details.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if float(line[5]) >= 0.5:
           over.append(line[5])
           count_over +=1
        else:
           under.append(line[5])
           count_under +=1    
