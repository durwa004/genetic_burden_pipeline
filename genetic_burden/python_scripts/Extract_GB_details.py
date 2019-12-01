import gzip
from statistics import mean  
from __future__ import division  

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
                
#Get unique GB  gene details
with open(path + "/unique_gb.txt", "r") as input_file, open(path + "/unique_gb_genes.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        print(line[9], file = output_file)

#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Create accession/gene symbol dictionary
genes = {}
with open(path + "/unique_gb.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        old = line[9].split(".")
        genes[old[0]] = "A"
with open(path + "/unique_gb_genes_with_symbols.txt", "r") as genes_file:
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       for key in genes.keys():
           if key == line1[0]:
               genes[key] = line1[1]

#Print out the unique variants and which horse they are present in
#0-12 not horses
with open(path + "/unique_gb.txt", "r") as input_file, open(path + "/unique_gb_individuals.txt", "w") as output_file:
    header = input_file.readline()
    header = header.split("\t")
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ind = []
        for i in range(len(line)):
            if "0/0" in line[i] or "./." in line[i]:
                next
            else:
                for key in horse_breed.keys():
                    if horse_breed[key] == "A":
                        next
                    elif header[i] == key:
                        a = key + ":" + horse_breed[key]
                        ind.append(a)
        c = line[9].split(".")
        print("\t".join(line[0:9]),genes[c[0]], line[10],
                  line[11], line[12], "\t".join(ind), sep = "\t", file = output_file)

#               
#Do for all GB variants
genes = {}
lof_genes = {}
with open(path + "/genetic_burden_details.txt", "r") as input_file, open(path + "/gb_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        genes[line[8]] = "A"
        print(line[8], file = output_file)
        if line[11] == "y":
            lof_genes[line[8]] = "A"
print(len(set(genes.keys())))
print(len(set(lof_genes.keys())))
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Create accession/gene symbol dictionary
with open(path + "/gb_genes_symbols.txt", "r") as genes_file:
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       for key in genes.keys():
           a = key.split(".")
           if a[0] == line1[0]:
               if line1[1] == "-":
                    genes[key] = key
               else:
                   genes[key] = line1[1]
       for key in lof_genes.keys():
           a = key.split(".")
           if a[0] == line1[0]:
               if line1[1] == "-":
                    lof_genes[key] = key
               else:
                   lof_genes[key] = line1[1]

print(len(set(genes.values())))
print(len(set(lof_genes.values())))

#Check if any missing
for key in genes.keys():
    if genes[key] == "A":
        print(key)
        
# number of genes affected by lof and GB variants
#Get number of variants per transcript
# number of lof variants per transcript
count = 0
AC = 0
AF = 0
lof_count = 0
lof_AC = 0
lof_AF = 0
with open(path + "/genetic_burden_details.txt", "r") as input_file:
    header = input_file.readline()
    header = header.split("\t")
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        count +=1
        AC += int(line[4])
        AF += float(line[5]) 
        if "," in genes[line[8]]:
            a = genes[line[8]] + "," + line[4]
            genes[line[8]] = a
        else:
            a = genes[line[8]]
            genes[line[8]] = a + "*" + line[4]            
        if line[11] == "y":
            lof_count +=1
            lof_AC += int(line[4])
            lof_AF += float(line[5])
            if "," in lof_genes[line[8]]:
                a = lof_genes[line[8]] + "," + line[4]
                lof_genes[line[8]] = a
            else:
                a = lof_genes[line[8]]
                lof_genes[line[8]] = a + "*" + line[4]              
                    
print(AC/count)
print(AF/count)
print(lof_AC/lof_count)
print(lof_AF/lof_count)


count_single = 0
count_multiple = 0
length = 0
n_variants = []
single_gene = {}
multiple_gene = {}
with open(path + "/multiple_variants_genes.txt", "w") as output_file:
    for key, value in genes.items():
        a = value.split(",")
        b = value.split("*")
        n_variants.append(len(a))
        if len(a) <2:
            count_single +=1
            if genes[key] in single_gene.keys():
                if genes[key] in multiple_gene.keys():
                    b = multiple_gene[genes[key]] + ":" + "1"
                    multiple_gene[genes[key]] = b
                else:
                    multiple_gene[genes[key]] = "1"
            else:
                single_gene[genes[key]] = "1"
        else:
            count_multiple +=1
            if genes[key] in multiple_gene.keys():
                b = multiple_gene[genes[key]] + ":" + "1"
                multiple_gene[genes[key]] = b
            else:
                multiple_gene[genes[key]] = "1"
            if len(a) > length:
                length = len(a)
            if len(a) > 5:
                print(genes[key], file = output_file)
            
print(len(single_gene))
print(len(multiple_gene))

len(n_variants) #5,821  
mean(n_variants) #1.52
n_variants_sort = sorted(n_variants)
print(n_variants_sort[0]) #1
print(n_variants_sort[-1]) #26
count_single
count_multiple

with open(path + "/gb_genes_symbols.txt", "r") as genes_file:
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       for key in genes.keys():
           if key == line1[0]:
               if line1[1] == "-":
                    genes[key] = key
               else:
                   genes[key] = line1[1]
       for key in lof_genes.keys():
           if key == line1[0]:
               if line1[1] == "-":
                    lof_genes[key] = key
               else:
                   lof_genes[key] = line1[1]
                   
len(genes)
len(set(genes.values()))
len(lof_genes)
len(set(lof_genes.values()))
#Do for LOF variants
