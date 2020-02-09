import os

#Get protein names from gene description details
unch = 0
desc = []
genes = {}
old = []
with open("lof_high_AF_gene_desc.txt", "r") as input_file, open("lof_high_AF_for_Uniprot.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[1].split("Gene Symbol: ")
        if "LOC1" in a[1]:
            b = line[1].split("[Description: ")
            c = b[1].split("]")
            d = c[0].split()
            if d[0] == "uncharacterized":
                unch +=1
            elif d[0] == "probable":
                genes[(d[-1])] = d[1]
            else:
                desc.append(c[0])
                old.append(line[0])
        else:
            genes[line[0]] = "NA"
    for i in range(len(desc)):
        print("/Users/durwa004/edirect/esearch -db gene -query ", '"', desc[i], ' AND Homo sapiens [ORGN]" | /Users/durwa004/edirect/esummary | grep -w "Name" > ', old[i],"_names.txt", sep = "", file = output_file)


#get gene names from the esearch output

ncbi_names = {}
for filename in os.listdir("../lof/"):
    if filename.endswith("_names.txt"):
        with open(filename, "r") as input_file:
            orig = filename.split("_names.txt")
            line = input_file.readline()
            line = line.rstrip("\n").split(">")
            if len(line) > 1:
                a = line[1].split("<")
                if a[0] in ncbi_names.keys():
                    ab = ncbi_names[a[0]] + ":" + orig[0]
                    ncbi_names[a[0]] = ab
                else:
                    ncbi_names[a[0]] = orig[0]

genes.update(ncbi_names)
with open("lof_high_AF_for_enrichment.txt", "w") as output_file:
    for key in genes.keys():
        print(key, file = output_file)


mart = []
gene = []
with open("mart_export_high_AF.txt", "r") as input_file, open("mart_high_AF_tidy.txt", "w") as output_file:
    print("Gene\tGO_term", file = output_file)
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if len(line) >8:
            a = line[8].split()
            for i in range(len(a)):
                if "olfact" in a[i]:
                    a = "olfactory_reception"
                    mart.append(a)
                else:
                    mart.append(line[8])
            gene.append(line[9])
            print(line[9], "_".join(a), sep = "\t", file = output_file)

genes_orig = []
with open("lof_high_AF_genes.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        genes_orig.append(line[0])
          
        
