from __future__ import division  
import gzip
from statistics import mean  

#This is a very long script and likely can be tidied up at some point!!
path = "/panfs/roc/groups/3/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof"

#Need to get header from vcf file
horse_breed = {}
with gzip.open(path + "/../joint_intersect_without_Prze/thesis_intersect.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[0] == "#CHROM":
            line = line[9:]
            for i in range(len(line)):
                horse_breed[line[i]] = "A"
            break

#Add in breed
with open(path + "/../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for key in horse_breed.keys():
            if line[0] == key:
                horse_breed[key] = line[1]
#0-8 are not horses

#N transcripts/genes affected
genes = {}
with open(path + "/lof_details_brief.txt", "r") as input_file, open(path + "/lof_genes.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[10] in genes.keys():
            pass
        else:
            genes[line[10]] = line[16]
            print(line[10], file = output_file)
        if line[16] in genes.keys():
            pass
        else:
            genes[line[16]] = line[10]
            print(line[16], file = output_file)

#Get number of transcripts affected
transcripts = list(genes.keys()) + list(genes.values())                
len(set(transcripts))

#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol


#Gene enrichment
genes = {}
with open(path + "/lof_genes_with_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genes.keys():
               genes[line1[0]][line1[0]] = "NA"
           else:
               genes[line1[0]] = {}
               genes[line1[0]][line1[0]] = "NA"
       else:
           if line1[1] in genes.keys():
               genes[line1[1]][line1[0]] = "NA"
           else:
               genes[line1[1]] = {}
               genes[line1[1]][line1[0]] = "NA"
print(len(genes))

count = 0
for key in genes.keys():
    for key1 in genes[key].keys():
        count +=1
print(count)

#Get information about the number of variants and AF of each variant
with open(path + "/lof_details_brief.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[10].split(".")
        for gene in genes.keys():
            for transcript in genes[gene].keys():
                if transcript == a[0]:
                    if "," in genes[gene][a[0]]:
                        b = genes[gene][a[0]] + line[6] + ":" + line[7] + ","
                        genes[gene][a[0]] = b
                    else:
                        b = line[6] + ":" + line[7] + "," 
                        genes[gene][a[0]] = b

#Get mean (range) number of variants per gene and the type of variants for each one
#Print out genes with over the average number of variants?
n_variants = []
single = 0
multiple = 0
lots = 0
AC_list = []
AF_list = []
v_type_count = []
with open(path + "/lof_multi-variant_genes.txt", "w") as output_file, open(path + "/lof_high_AF_genes.txt", "w") as high_AF:
    for key,value in genes.items():
        variants = 0
        AF = []
        v_type = []
        AF_n = 0
        transcripts = []
        count = 0
        for item,value1 in genes[key].items():
            if value1 != "NA":
                a = value1.split(",")
                for i in range(len(a)):
                    if len(a[i]) > 1:
                        count +=1
                        transcripts.append(item)
                        ab = a[i].split(":")
                        if len(ab) >1:
                            variants +=1
                            v_type.append(ab[1])
                            AF_n += float(ab[0])
                            AF.append(ab[0])
        if count > 0:
            if count == 1:
                single +=1
            elif count >1:
                multiple +=1
            if count > 5:
                lots +=1
                print(key, file = output_file)
            SNP = v_type.count("SNP")
            indel = v_type.count("indel")
            S_i = str(SNP) + ":" + str(indel)
            v_type_count.append(S_i)
            n_variants.append(variants)
            ef = float(AF_n)/int(variants)
            AF_list.append(ef)
            if float(ef) > 0.50:
                print(key, file = high_AF)

print(single)
print(multiple)
print(lots)

n_variants_sort = sorted(n_variants)
print(n_variants_sort[0])
print(mean(n_variants_sort))
print(n_variants_sort[-1])

AF_list_sort = sorted(AF_list)
print(mean(AF_list_sort))
print(AF_list_sort[0])
print(AF_list_sort[-1])



#############################################################################
#Create accession/gene symbol dictionary for homozygous/heterozygous variants
genotype = {}
with open(path + "/lof_genes_with_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genotype.keys():
               genotype[line1[0]][line1[0]] = "NA"
           else:
               genotype[line1[0]] = {}
               genotype[line1[0]][line1[0]] = "NA"
       else:
           if line1[1] in genotype.keys():
               genotype[line1[1]][line1[0]] = "NA"
           else:
               genotype[line1[1]] = {}
               genotype[line1[1]][line1[0]] = "NA"
print(len(genotype))

count = 0
for key in genotype.keys():
    for key1 in genotype[key].keys():
        count +=1
print(count)

#Get hom/het details
with open(path + "/lof_details.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[10].split(".")
        for gene in genotype.keys():
            for transcript in genotype[gene].keys():
                count_hom = 0
                count_het = 0
                count_miss = 0
                count_ref = 0
                if transcript == a[0]:
                    genotypes = line[19:]
                    if "," in genotype[gene][a[0]]:
                        for item,value in enumerate(genotypes):
                            if "ref" == value:
                                count_ref +=1
                            elif "hom" == value:
                                count_hom +=1
                            elif "het" == value:
                                count_het +=1
                            elif "missing" == value:
                                count_miss +=1
                            else:
                                print(value)
                        b = genotype[gene][a[0]] + str(count_ref) + ":" + str(count_hom) + ":" + str(count_het) + ":" + str(count_miss) + ","
                        genotype[gene][a[0]] = b
                    else:
                        for item,value in enumerate(genotypes):
                            if "ref" == value:
                                count_ref +=1
                            elif "hom" == value:
                                count_hom +=1
                            elif "het" == value:
                                count_het +=1
                            elif "missing" == value:
                                count_miss +=1
                            else:
                                print(value)
                        b = str(count_ref) + ":" + str(count_hom) + ":" + str(count_het) + ":" + str(count_miss) + ","
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
with open(path + "/multiple_homozygous_lof_genotypes", "w") as output2:
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
        n_hom.append(trans_hom)
        lof_n_hom.append(trans_lof_hom)
        if count_hom > 0:
            transcript_hom.append(count_hom)
        if count_hom >5:
            lots.append(count_hom)
            print(key, file = output2)

    
print(len(transcript_hom))
print(len(lots))

n_hom_sort = sorted(n_hom)
print(mean(n_hom))
print(n_hom_sort[0])
print(n_hom_sort[-1])

transcript_hom_sort = sorted(transcript_hom)
print(mean(transcript_hom))
print(transcript_hom_sort[0])
print(transcript_hom_sort[-1])
                    
