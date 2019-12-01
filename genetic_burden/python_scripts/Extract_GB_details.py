import gzip
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
transcripts = []
lof_transcripts = []
with open(path + "/genetic_burden_details.txt", "r") as input_file, open(path + "/gb_genes.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        transcripts.append(line[8])
        print(line[8], file = output_file)
        if line[11] == "y":
            lof_transcripts.append(line[8])
print(len(set(transcripts)))
print(len(set(lof_transcripts)))
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Create accession/gene symbol dictionary
genes = {}
lof_genes = {}
with open(path + "/genetic_burden_details.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        old = line[8].split(".")
        genes[old[0]] = "A"
        if line[11] == "y":
            lof_genes[old[0]] = "A"
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


#Get number of variants per gene
# number of lof variants per gene
# number of genes affected by lof variants
count = 0
AC = 0
AF = 0
lof_count = 0
lof_AC = 0
lof_AF = 0

AC_dict = {}
AF_dict = {}
with open(path + "/genetic_burden_details.txt", "r") as input_file:
    header = input_file.readline()
    header = header.split("\t")
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        count +=1
        AC += int(line[4])
        AF += float(line[5])
        if line[11] == "y":
            lof_count +=1
            lof_AC += int(line[4])
            lof_AF += float(line[5])
        for key,value in genes.items():
            if key in line[8]:
                if key in AC_dict.keys():
                    a = AC_dict[key]
                    c = int(a) + int(line[4])
                    ac = c/2
                    b = AF_dict[key]
                    d = float(b) + float(line[5])
                    bd = d/2
                    AC_dict[key] = ac
                    AF_dict[key] = bd
                else:
                    AC_dict[key] = str(line[4])
                    AF_dict[key] = str(line[5])
                    
print(AC/count)
print(AF/count)
print(lof_AC/lof_count)
print(lof_AF/lof_count)

a = 0
for key, value in AC_dict.items():
    if int(value) < a:
        a = value

#Do for LOF variants
