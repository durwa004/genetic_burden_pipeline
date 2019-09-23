import gzip
import os

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open("SnpEff_coding.vcf.gz", "rt") as input_file, open("SnpEff_coding_tidy.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line:
            print("\t".join(line[0:7]), "AC", "AF", "Consequence","Impact", "Gene", "Where", "coding_change", "protein_change", "\t".join(line[9:]), file = output_file, sep = "\t")
        elif "#" in line[0]:
            next
        else:
            ab = line[7].split("ANN=")
            bc = ab[1].split("|")
            if bc[0] == line[4]:
                a = line[7].split("AC=")
                b = a[1].split(";")
                AC = b[0]
                if "AF_" in line[7]:
                    c = line[7].split("AF=")
                    d = c[1].split(";")
                    AF = d[0]
                else:
                    AF = "NA"
                consequence = bc[1]
                impact = bc[2]
                gene = bc[3]
                where = bc[7]
                coding = bc[9]
                protein = bc[10]
            else:
                a = line[7].split("AC=")
                b = a[1].split(";")
                cd = b[0].split(",")
                AC =cd[0]
                if "AF_" in line[7]:
                    c = line[7].split("AF=")
                    d = c[1].split(";")
                    de = d[0].split(",")
                    AF = de[0]
                else:
                    AF = "NA"
                consequence = bc[1]
                impact = bc[2]
                gene = bc[3]
                where = bc[7]
                coding = bc[9]
                protein = bc[10]
            print("\t".join(line[0:7]), AC, AF, consequence, impact, gene, where, coding, protein, "\t".join(line[9:]), file = output_file, sep = "\t")

#Extract annovar variants into same format as snpeff variants
with open("annovar_coding.txt", "r") as input_file, open("SnpEff_coding_tidy.txt", "r") as header, open("annovar_coding_tidy.txt", "w") as output_file:
    info = header.readline()
    print(info, file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        consequence = line[1]
        if "frameshift" in line[1] or "stopgain" in line[1] or "stoplost" in line[1]:
            impact = "HIGH"
        elif "nonframeshift" in line[1] or "nonsynonymous" in line[1]:
            impact = "MODERATE"
        elif "synonymous SNV" in line[1]:
            impact = "LOW"
        elif "unknown" in line[1]:
            impact = "UNKNOWN"
        if line[2] == "UNKNOWN":
            b = line[18].split("AC=")
            c = b[1].split(";")
            AC = c[0]
            e = line[18].split("AF=")
            f = e[1].split(";")
            AF = f[0]
            print("\t".join(line[11:18]), AC, AF, consequence, impact, "NA", "NA", "NA", "NA", "\t".join(line[20:]), file = output_file, sep = "\t")
        else:
            a = line[2].split(":")
            gene = a[0]
            coding = a[3]
            protein = a[4]
            b = line[18].split("AC=")
            c = b[1].split(";")
            AC = c[0]
            e = line[18].split("AF=")
            f = e[1].split(";")
            AF = f[0]
            print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein, "\t".join(line[20:]), file = output_file, sep = "\t")
 

         
#Need to split annovar/snpeff union file by breed for analysis:
horse = []
breed = []
with open("../ibio_horses_with_breeds.txt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse.append(line[0])      
        if line[1] == "Arabian":
            breed.append("Arabian")
        elif line[1] == "Belgian":
            breed.append("Belgian")
        elif line[1] == "Clydesdale":
            breed.append("Clydesdale")
        elif line[1] == "Icelandic" or line[1] == "Icelandic Horse":
            breed.append("Icelandic")
        elif line[1] == "Morgan":
            breed.append("Morgan")
        elif line[1] == "QH" or line[1] == "QH (App)" or line[1] == "Quarter Horse":
            breed.append("QH")
        elif line[1] == "Shetland":
            breed.append("Shetland")
        elif line[1] =="STB" or "StandardBred" in line[1]:
            breed.append("STB")
        elif line[1] == "Thoroughbred":
            breed.append("TB")
        elif line[1] == "Welsh Pony":
            breed.append("WP")
        else:
            breed.append("Other")
#Print out tidied up breed list:
with open("ibio_ids_with_breed_groups.txt", "w") as output_file:
    for i in range(len(breed)):
        print(horse[i], breed[i], file = output_file,sep="\t")

#Print out ids of different breeds
for i,v in enumerate(set(breed)):
    with open(v + "_ids.txt", "w") as output_file:
        for item in range(len(breed)):
           if breed[item] == v:
               print(horse[item],file=output_file)
           else:
               continue

data = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect"   
with gzip.open(data + "/SnpEff/SnpEff_annovar_coding_union.vcf.gz", "rt") as input_file, open(data + "/snpeff_annovar_union_summary.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "##" in line[0]:
            input_file.readline()
        elif "#CHROM" in line[0]:
            print("CHROM\tPOS\tLOF\tAF\tGene\t","\t".join(line[9:]), file = output_file,sep="")
            header = ["CHROM","POS","LOF","AF","Gene"] + line[9:]
            with open("ibio_ids_with_breed_groups.txt", "r") as breed_file:
                breed_line = {}
                breed = []
                for line in breed_file:
                    line = line.rstrip("\n").split("\t")
                    breed_line[line[0]] = line[1]
                    breed_line["CHROM"] = "NA"
                    breed_line["POS"] = "NA"
                    breed_line["LOF"] = "NA"
                    breed_line["AF"] = "NA"
                    breed_line["Gene"] = "NA"
                    breed_line["TWILIGHT"] = "TB"
            for i in range(len(header)):
                a = breed_line[header[i]]
                breed.append(a)
            print("\t".join(breed),file = output_file)
        else:
            if "LOF" in line[7]:
                lof = "lof"
            else:
                lof = "no"
            info = line[7].split(";")
            if "AF=" in info[1]:
                if "," in info[1]:
                    b = info[1].split(",")
                    AF = b[0].split("AF=")
                    AF = AF[1]
                else:
                    AF = info[1].split("AF=")
	            AF = AF[1]
            else:
                AF = info[1]
            ANN = line[7].split("ANN")
            ANN = ANN[1]
            a = ANN.split("|")
            effect = a[1]
            if "gene-" in line[7]:
                c = line[7].split("gene-")
                gene = c[1].split("|")
                gene = gene[0]
            else: 
                gene = a[3]
            IDs = line[9:]
            for i,v in enumerate(IDs):
                if "1/1" in v:
                    IDs[i] = "2"
                elif "0/1" in v:
                    IDs[i] = "1"
                elif "0/0" in v:
                    IDs[i] = "0"
                else:
                    IDs[i] = "NA"
            print(line[0],line[1], lof, AF,gene,"\t".join(IDs),file = output_file,sep = "\t")



#Having problems with the pandas df = try without so I can get some results!
with open("snpeff_annovar_union_summary.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")

#Read file back in and figure out statistics
import statistics
import numpy as np
import pandas as pd

high_df = pd.read_table("snpeff_annovar_union_summary.txt", sep="\t")

#Cannot figure this out. So write to my laptop and use R if necessary
lof = high_df[high_df["LOF"] == "lof"]

f
#Transpose table
no_high = high_df_T["CHROM"]["0/1"].value_counts()
high_df_T = high_df.T

lof = high_df[(high_df["LOF"] == "lof") | (high_df["LOF"] != "no")]
lof.to_csv("ibio_lof_union.csv")






#Number of homozgyous high impact variants
high_hom = high_df[(high_df[6:] == "1/1")
data_hom = data[(data["1KG16"] == "1/1") & (data["1KG17"] == "1/1") & (data["1KG18"] == "1/1")]

mean_AF = high_df["AF"].value_counts()

no_high = high_df["ID"].value_counts()
statistics.mean(no_high) #2343
statistics.median(no_high)#2304
min(no_high)#348
max(no_high)#3468

#Number of genes affected
high_gene_totals = high_df["ANN_gene"].value_counts()
len(high_gene_totals) #8746
statistics.mean(high_gene_totals) #61
statistics.median(high_gene_totals)#3
min(high_gene_totals)#1
max(high_gene_totals)#2028
#Get list of genes
total_genes = list(set(high_df["ANN_gene"]))


#Number of LOF variants per individual:
LOF = high_df["LOF?"] == "y"  
LOF_no = high_df[LOF]  
LOF_totals = LOF_no["ID"].value_counts()
statistics.mean(LOF_totals) #1888
statistics.median(LOF_totals) #1863
min(LOF_totals) #246
max(LOF_totals) #2733

#No. of genes:
LOF_gene_totals = LOF_no["ANN_gene"].value_counts()
len(LOF_gene_totals) #7049
statistics.mean(LOF_gene_totals) #61
statistics.median(LOF_gene_totals)#2
min(LOF_gene_totals)#1
max(LOF_gene_totals)#2026

#Calculate average number of homozygous high impact mutations per individual
hom = high_df["het/hom"] == "1/1"
hom_no = high_df[hom]
hom_totals = hom_no["ID"].value_counts()
statistics.mean(hom_totals) #1952
statistics.median(hom_totals) #1909
min(hom_totals) #268
max(hom_totals)#2825

#Calculate average number of homozygous LOF mutations per individual
LOF_hom = LOF_no["het/hom"] == "1/1"
LOF_hom_n = LOF_no[LOF_hom]
LOF_hom_totals = LOF_hom_n["ID"].value_counts()
statistics.mean(LOF_hom_totals) #1626
statistics.median(LOF_hom_totals) #1597
min(LOF_hom_totals) #194
max(LOF_hom_totals) #2331
#Get list of genes
LOF_total_genes = list(set(LOF_no["ANN_gene"]))

#Number of genes affected:
hom_gene_totals = hom_no["ANN_gene"].value_counts()
len(hom_gene_totals) #3874
statistics.mean(hom_gene_totals) #115
statistics.median(hom_gene_totals)#73
min(hom_gene_totals)#1
max(hom_gene_totals)#1734
#Get list of genes:
hom_genes = list(set(hom_no["ANN_gene"]))


#Number of genes affected by LOF
LOF_gene_totals = LOF_hom_n["ANN_gene"].value_counts()
len(LOF_gene_totals) #2944
statistics.mean(LOF_gene_totals) #127
statistics.median(LOF_gene_totals)#103
min(LOF_gene_totals)#1
max(LOF_gene_totals)#1730

#Get list of genes:
LOF_hom_genes = list(set(LOF_hom_n["ANN_gene"]))

#Print high impact and LOF genes to new files for gene ontology:
#high impact all
with open("230_horses_high_impact_genes.txt", "w") as output_file:
    for i in range(len(total_genes)):
        print(total_genes[i], file = output_file, sep = "\n")

#LOF all
with open("230_horses_LOF_genes.txt", "w") as output_file:
    for i in range(len(LOF_total_genes)):
        print(LOF_total_genes[i], file = output_file, sep = "\n")

#high impact hom
with open("230_horses_hom_high_impact_genes.txt", "w") as output_file:
    for i in range(len(hom_genes)):
        print(hom_genes[i], file = output_file, sep = "\n")

#LOF hom
with open("230_horses_hom_LOF_genes.txt", "w") as output_file:
    for i in range(len(LOF_hom_genes)):
        print(LOF_hom_genes[i], file = output_file, sep = "\n")


#Need to add in breed as well.
#High impact variants by breed:
all_high = high_df["Breed_group"].value_counts()

#Determine if there is a statistical difference between breeds
import scipy
from scipy import stats
no_high = high_df["ID"].value_counts()
high_dict = dict(no_high)

#Take ID_breed and breed_group from earlier
#All high impact variants
high_ID = list()
high_no = list()
high_breed = list()

for key, value in enumerate(high_dict.keys()):
    for i in range(len(ID_breed)):
        if ID_breed[i] == value:
            high_ID.append(ID_breed[i])
            high_no.append(high_dict.get(value))
            high_breed.append(breed_group[i])
            
stats.kruskal(high_no, high_breed)
#p = 4.129 x 10e-7

#Double check these results in R:
with open("230_horses_high_impact_by_breed.txt", "w") as output_file:
    for i in range(len(high_ID)):
        print(high_ID[i], high_no[i], high_breed[i], sep = "\t", file = output_file)

#Just homozygous high impact variants
high_hom_dict = dict(hom_totals)

high_hom_ID = list()
high_hom_no = list()
high_hom_breed = list()

for key, value in enumerate(high_hom_dict.keys()):
    for i in range(len(ID_breed)):
        if ID_breed[i] == value:
            high_hom_ID.append(ID_breed[i])
            high_hom_no.append(high_hom_dict.get(value))
            high_hom_breed.append(breed_group[i])
            
stats.kruskal(high_hom_no, high_hom_breed)
#p = 4.129 x 10e-77

#Double check these results in R:
with open("230_horses_hom_high_impact_by_breed.txt", "w") as output_file:
    for i in range(len(high_hom_ID)):
        print(high_hom_ID[i], high_hom_no[i], high_hom_breed[i], sep = "\t", file = output_file)

#All LOF variants

LOF_dict = dict(LOF_totals)
LOF_ID = list()
LOF_no = list()
LOF_breed = list()

for key, value in enumerate(LOF_dict.keys()):
    for i in range(len(ID_breed)):
        if ID_breed[i] == value:
            LOF_ID.append(ID_breed[i])
            LOF_no.append(LOF_dict.get(value))
            LOF_breed.append(breed_group[i])
            
#Double check these results in R:
with open("230_horses_LOF_by_breed.txt", "w") as output_file:
    for i in range(len(LOF_ID)):
        print(LOF_ID[i], LOF_no[i], LOF_breed[i], sep = "\t", file = output_file)

#Just homozygous LOF variants
LOF_hom_dict = dict(LOF_hom_totals)

LOF_hom_ID = list()
LOF_hom_no = list()
LOF_hom_breed = list()

for key, value in enumerate(LOF_hom_dict.keys()):
    for i in range(len(ID_breed)):
        if ID_breed[i] == value:
            LOF_hom_ID.append(ID_breed[i])
            LOF_hom_no.append(LOF_hom_dict.get(value))
            LOF_hom_breed.append(breed_group[i])
            

#Double check these results in R:
with open("230_horses_hom_LOF_by_breed.txt", "w") as output_file:
    for i in range(len(LOF_hom_ID)):
        print(LOF_hom_ID[i], LOF_hom_no[i], LOF_hom_breed[i], sep = "\t", file = output_file)


#Want to compare the number of high impact variants by breed
#Need 2 columns: Number of high impact variants per individual and then the breed for each individual:
#Number of high impact variants per individual:
no_high = high_df["ID"].value_counts()
high_per_ind = list(no_high)

#QH
QH = high_df["Breed_group"] == "QH"
QH_df = high_df[QH]

QH_high = QH_df["ID"].value_counts()
statistics.mean(QH_high)
statistics.median(QH_high)
min(QH_high)
max(QH_high)

#Arab
Arab = high_df["Breed_group"] == "Arab"
Arab_df = high_df[Arab]

Arab_high = Arab_df["ID"].value_counts()
statistics.mean(Arab_high)
statistics.median(Arab_high)
min(Arab_high)
max(Arab_high)

#Belgian
Belgian = high_df["Breed_group"] == "Belgian"
Belgian_df = high_df[Belgian]

Belgian_high = Belgian_df["ID"].value_counts()
statistics.mean(Belgian_high)
statistics.median(Belgian_high)
min(Belgian_high)
max(Belgian_high)

#Clydesdale
Clydesdale = high_df["Breed_group"] == "Clydesdale"
Clydesdale_df = high_df[Clydesdale]

Clydesdale_high = Clydesdale_df["ID"].value_counts()
statistics.mean(Clydesdale_high)
statistics.median(Clydesdale_high)
min(Clydesdale_high)
max(Clydesdale_high)

#Morgan
Morgan = high_df["Breed_group"] == "Morgan"
Morgan_df = high_df[Morgan]

Morgan_high = Morgan_df["ID"].value_counts()
statistics.mean(Morgan_high)
statistics.median(Morgan_high)
min(Morgan_high)
max(Morgan_high)

#STB
STB = high_df["Breed_group"] == "STB"
STB_df = high_df[STB]

STB_high = STB_df["ID"].value_counts()
statistics.mean(STB_high)
statistics.median(STB_high)
min(STB_high)
max(STB_high)

#TB
TB = high_df["Breed_group"] == "TB"
TB_df = high_df[TB]

TB_high = TB_df["ID"].value_counts()
statistics.mean(TB_high)
statistics.median(TB_high)
min(TB_high)
max(TB_high)

#WP
WP = high_df["Breed_group"] == "WP"
WP_df = high_df[WP]

WP_high = WP_df["ID"].value_counts()
statistics.mean(WP_high)
statistics.median(WP_high)
min(WP_high)
max(WP_high)

#Other
Other = high_df["Breed_group"] == "Other"
Other_df = high_df[Other]

Other_high = Other_df["ID"].value_counts()
statistics.mean(Other_high)
statistics.median(Other_high)
min(Other_high)
max(Other_high)
np.std(Other_high)
len(Other_high)
