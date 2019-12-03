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
                
#Get unique GB  gene details
with open(path + "/unique_gb.txt", "r") as input_file, open(path + "/unique_gb_genes.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        print(line[9], file = output_file)

#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Create accession/gene symbol dictionary
genes_unique = {}
with open(path + "/unique_gb.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        old = line[9].split(".")
        genes_unique[old[0]] = "A"
with open(path + "/unique_gb_genes_with_symbols.txt", "r") as genes_file:
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       for key in genes_unique.keys():
           if key == line1[0]:
               genes_unique[key] = line1[1]

#Print out the unique variants and which horse they are present in
#0-12 not horses
with open(path + "/unique_gb.txt", "r") as input_file, open(path + "/unique_gb_individuals.txt", 
         "w") as output_file, open(path + "/unique_gb_genes.txt", "w") as output2:
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
        print("\t".join(line[0:9]),genes_unique[c[0]], line[10],
                  line[11], line[12], "\t".join(ind), sep = "\t", file = output_file)
        print(genes_unique[c[0]], file = output2)

#               
#Do for all GB variants
with open(path + "/genetic_burden_details.txt", "r") as input_file, open(path + "/gb_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        print(line[8], file = output_file)
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Create accession/gene symbol dictionary
genes = {}
with open(path + "/gb_genes_symbols.txt", "r") as genes_file:
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

#Add in AC and AF to gene/transcript dictionary
with open(path + "/genetic_burden_details.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[8].split(".")
        for gene in genes.keys():
            for transcript in genes[gene].keys():
                if transcript == a[0]:
                    if "," in genes[gene][a[0]]:
                        b = genes[gene][a[0]] + line[0] + ":" + line[1] + ":" + line[4] + ":" + line[5] + ":" + line[11] + ","
                        genes[gene][a[0]] = b
                    else:
                        b = line[0] + ":" + line[1] + ":" + line[4] + ":" + line[5] + ":" + line[11] + ","
                        genes[gene][a[0]] = b
                        
# number of genes affected by GB variants
n_variants = []
single = 0
multiple = 0
lots = 0
AC_list = []
AF_list = []
with open(path + "/exploring_GB.txt", "w") as output_file, open(path + "/multiple_variants_genes.txt", "w") as output2:
    print("gene\ttranscript\tlocation\tAC\tAF\tLOF", file = output_file)
    for key,value in genes.items():
        variants = 0
        location = []
        AC_n = 0
        AF_n = 0
        AC = []
        AF = []
        LOF = []
        transcripts = []
        count = 0
        for item,value1 in genes[key].items():
            a = value1.split(",")
            count +=1
            for i in range(len(a)):
                if a[i] == "":
                    next
                else:
                    count +=1
                    transcripts.append(item)
                    ab = a[i].split(":")
                    b = ab[0] + ":" +ab[1]
                    variants +=1
                    location.append(b)
                    AC_n += int(ab[2])
                    AF_n += float(ab[3])
                    AC.append(ab[2])
                    AF.append(ab[3])
                    LOF.append(ab[4])
        if (count-1) == 1:
            single +=1
        else:
            multiple +=1
        if (count-1) > 5:
            lots +=1
            print(key, file = output2)
        n_variants.append(variants)
        cd = int(AC_n)/len(location)
        AC_list.append(cd)
        ef = float(AF_n)/len(location)
        AF_list.append(ef)
        for i in range(len(location)):
            print(key, transcripts[i], location[i], AC[i], AF[i], LOF[i], sep = "\t", file = output_file)

print(single)
print(multiple)
print(lots)

n_variants_sort = sorted(n_variants)
n_variants_sort[0]
n_variants_sort[-1]

AC_list_sort = sorted(AC_list)
AC_list_sort[0]
AC_list_sort[-1]

AF_list_sort = sorted(AF_list)
AF_list_sort[0]
AF_list_sort[-1]

#Get number of variants per transcript
n_variants = []
AC_list = []
AF_list = []
single = 0
multiple = 0
lots = 0
for key,value in genes.items():
    for item,value1 in genes[key].items():
        variants = 0
        AC_n = 0
        AF_n = 0
        a = value1.split(",")
        for i in range(len(a)):
            if a[i] == "":
                next
            else:
                ab = a[i].split(":")
                b = ab[0] + ":" +ab[1]
                variants +=1
                AC_n += int(ab[2])
                AF_n += float(ab[3])
            n_variants.append(variants)
            cd = int(AC_n)/(len(a)-1)
            AC_list.append(cd)
            ef = float(AF_n)/(len(a)-1)
            AF_list.append(ef)
        if len(a) <3:
            single +=1
        else:
            multiple +=1
        if len(a) >5:
            lots +=1
print(single)
print(multiple)
print(lots)

n_variants_sort = sorted(n_variants)
print(mean(n_variants))
print(n_variants_sort[0])
print(n_variants_sort[-1])

AC_list_sort = sorted(AC_list)
print(mean(AC_list))
print(AC_list_sort[0])
print(AC_list_sort[-1])

AF_list_sort = sorted(AF_list)
print(mean(AF_list))
print(AF_list_sort[0])
print(AF_list_sort[-1])

#############################################################################
# Get number of genes affected by lof 
# number of lof variants per transcript

n_variants = []
single = 0
multiple = 0
lots = 0
AC_list = []
AF_list = []
with open(path + "/exploring_LOF.txt", "w") as output_file, open(path + "/multiple_variants_lof_genes.txt", "w") as output2:
    print("gene\ttranscript\tlocation\tAC\tAF\tLOF", file = output_file)
    for key,value in genes.items():
        variants = 0
        location = []
        AC_n = 0
        AF_n = 0
        AC = []
        AF = []
        LOF = []
        transcripts = []
        count = 0
        for item,value1 in genes[key].items():
            a = value1.split(",")
            for i in range(len(a)):
                if a[i] == "":
                    next
                else:
                    ab = a[i].split(":")
                    if ab[4] == "y":
                        count +=1
                        transcripts.append(item)
                        b = ab[0] + ":" +ab[1]
                        variants +=1
                        location.append(b)
                        AC_n += int(ab[2])
                        AF_n += float(ab[3])
                        AC.append(ab[2])
                        AF.append(ab[3])
                        LOF.append(ab[4])
        if count >0:
            if (count) == 1:
                single +=1
            else:
                multiple +=1
            if (count) > 5:
                lots +=1
                print(key, file = output2)
            n_variants.append(variants)
            cd = int(AC_n)/len(location)
            AC_list.append(cd)
            ef = float(AF_n)/len(location)
            AF_list.append(ef)
            for i in range(len(location)):
                print(key, transcripts[i], location[i], AC[i], AF[i], LOF[i], sep = "\t", file = output_file)

print(single)
print(multiple)
print(lots)

n_variants_sort = sorted(n_variants)
print(mean(n_variants_sort))
print(n_variants_sort[0])
print(n_variants_sort[-1])

AC_list_sort = sorted(AC_list)
print(mean(AC_list))
print(AC_list_sort[0])
print(AC_list_sort[-1])

AF_list_sort = sorted(AF_list)
print(mean(AF_list))
print(AF_list_sort[0])
print(AF_list_sort[-1])

#Get number of variants per transcript
n_variants = []
AC_list = []
AF_list = []
single = 0
multiple = 0
lots = 0
for key,value in genes.items():
    count = 0
    for item,value1 in genes[key].items():
        variants = 0
        AC_n = 0
        AF_n = 0
        a = value1.split(",")
        for i in range(len(a)):
            if a[i] == "":
                next
            else:
                ab = a[i].split(":")
                if ab[4] == "y":
                    count +=1
                    transcripts.append(item)
                    b = ab[0] + ":" +ab[1]
                    variants +=1
                    location.append(b)
                    variants +=1
                    AC_n += int(ab[2])
                    AF_n += float(ab[3])
            n_variants.append(variants)
            cd = int(AC_n)/(len(a)-1)
            AC_list.append(cd)
            ef = float(AF_n)/(len(a)-1)
            AF_list.append(ef)
        if count >0:
            if count == 1:
                single +=1
            else:
                multiple +=1
            if count >5:
                lots +=1

print(single)
print(multiple)
print(lots)

n_variants_sort = sorted(n_variants)
print(mean(n_variants))
print(n_variants_sort[0])
print(n_variants_sort[-1])

AC_list_sort = sorted(AC_list)
print(mean(AC_list))
print(AC_list_sort[0])
print(AC_list_sort[-1])

AF_list_sort = sorted(AF_list)
print(mean(AF_list))
print(AF_list_sort[0])
print(AF_list_sort[-1])


#############################################################################
#Create accession/gene symbol dictionary for homozygous/heterozygous variants
genotype = {}
with open(path + "/gb_genes_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genotype.keys():
               genotype[line1[0]][line1[0]] = {}
           else:
               genotype[line1[0]] = {}
               genotype[line1[0]][line1[0]] = {}
       else:
           if line1[1] in genotype.keys():
               genotype[line1[1]][line1[0]] = {}
           else:
               genotype[line1[1]] = {}
               genotype[line1[1]][line1[0]] = {}
print(len(genotype))

count = 0
for key in genotype.keys():
    for key1 in genotype[key].keys():
        count +=1
print(count)

#Add in AC and AF to gene/transcript dictionary
with open(path + "/genetic_burden_details.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[8].split(".")
        for gene in genotype.keys():
            for transcript in genotype[gene].keys():
                count_hom = 0
                count_het = 0
                count_miss = 0
                count_ref = 0
                if transcript == a[0]:
                    genotypes = line[13:]
                    if "," in genotype[gene][a[0]]:
                        for item,value in enumerate(genotypes):
                            if "0/0" in value:
                                count_ref +=1
                            elif "1/1" in value or "1/2" in value or "1/3" in value or "2/2" in value or "3/3" in value or "2/1" in value or "2/3" in value:
                                count_hom +=1
                            elif "0/1" in value or "0/2" in value or "0/3" in value:
                                count_het +=1
                            elif "./." in value:
                                count_miss +=1
                            else:
                                print(value)
                        b = genotype[gene][a[0]] + str(count_ref) + ":" + str(count_hom) + ":" + str(count_het) + ":" + str(count_miss) + ":" + line[11] + ","
                        genotype[gene][a[0]] = b
                    else:
                        for item,value in enumerate(genotypes):
                            if "0/0" in value:
                                count_ref +=1
                            elif "1/1" in value or "1/2" in value or "1/3" in value or "2/2" in value or "3/3" in value or "2/1" in value or "2/3" in value:
                                count_hom +=1
                            elif "0/1" in value or "0/2" in value or "0/3" in value:
                                count_het +=1
                            elif "./." in value:
                                count_miss +=1
                            else:
                                print(value)
                        b = str(count_ref) + ":" + str(count_hom) + ":" + str(count_het) + ":" + str(count_miss) + ":" + line[11] + ","
                        genotype[gene][a[0]] = b

#Get number of homozygotes/heterozygotes for each gene
n_hom = []
n_het = []
lof_n_hom = []
lof_n_het = []
transcript_hom = []
transcript_het = []
lots = []
lof_transcript_hom = []
lof_transcript_het = []
lof_lots = []
with open(path + "/multiple_homozygous_gb_genotypes.txt", "w") as output_file, open(path+"/multiple_homozygous_lof_genotypes", "w") as output2:
    for key,value in genotype.items():
        count_hom = 0
        count_het = 0
        lof_hom = 0
        lof_het = 0
        trans_hom = 0
        trans_lof_hom = 0
        for item,value1 in genotype[key].items():
            a = value1.split(",")
            for i in range(len(a)):
                ab = a[i].split(":")
                if len(ab) >1:
                    if int(ab[1]) >0:
                        count_hom += 1
                        trans_hom += int(ab[1]) 
                    if ab[4] == "y":
                        if int(ab[1]) > 0:
                            lof_hom += 1
                            trans_lof_hom += int(ab[1])
        n_hom.append(trans_hom)
        lof_n_hom.append(trans_lof_hom)
        if count_hom > 0:
            transcript_hom.append(count_hom)
        if count_hom >5:
            lots.append(count_hom)
            print(key, file = output_file)
        if lof_hom > 0:
            lof_transcript_hom.append(count_hom)
        if lof_hom >5:
            lof_lots.append(count_hom)
            print(key, file = output2)

    
print(len(transcript_hom))
print(len(lots))

n_hom_sort = sorted(n_hom)
print(mean(n_hom))
print(n_hom_sort[0])
print(n_hom_sort[-1])

print(len(lof_transcript_hom))
print(len(lof_transcript_het))
print(len(lof_lots))

lof_hom_sort = sorted(lof_n_hom)
print(mean(lof_hom_sort))
print(lof_hom_sort[0])
print(lof_hom_sort[-1])
                    
#Get number of homozygotes/heterozygotes for each transcript
n_hom = []
n_het = []
lof_n_hom = []
lof_n_het = []
transcript_hom = []
transcript_het = []
lots = []
lof_transcript_hom = []
lof_transcript_het = []
lof_lots = []
for key,value in genotype.items():
    count = 0
    for item,value1 in genotype[key].items():
        count_hom = 0
        count_het = 0
        lof_hom = 0
        lof_het = 0
        a = value1.split(",")
        for i in range(len(a)):
            ab = a[i].split(":")
            if len(ab) >1:
                if int(ab[1]) >0:
                    count_hom += 1
                    n_hom.append(int(ab[1]))
                elif int(ab[2]) > 0:
                    count_het +=1
                    n_het.append(ab[1])  
                if ab[4] == "y":
                    if int(ab[1]) > 0:
                        lof_hom += 1
                        lof_n_hom.append(int(ab[1]))
                    elif int(ab[2]) > 0:
                        lof_het +=1
                        lof_n_het.append(ab[2])
        if count_hom > 0:
            transcript_hom.append(count_hom)
        if count_hom >5:
            lots.append(count_hom)
        if count_het >0:
            transcript_het.append(count_het)
        if lof_hom > 0:
            lof_transcript_hom.append(count_hom)
        if lof_hom >5:
            lof_lots.append(count_hom)
        if lof_het >0:
            lof_transcript_het.append(count_het) 

    
print(len(transcript_hom))
print(len(transcript_het))
print(len(lots))

n_hom_sort = sorted(n_hom)
print(mean(n_hom))
print(n_hom_sort[0])
print(n_hom_sort[-1])

print(len(lof_transcript_hom))
print(len(lof_transcript_het))
print(len(lof_lots))

lof_hom_sort = sorted(lof_n_hom)
print(mean(lof_hom_sort))
print(lof_hom_sort[0])
print(lof_hom_sort[-1])
                    