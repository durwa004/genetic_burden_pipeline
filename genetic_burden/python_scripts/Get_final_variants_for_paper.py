path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper/Final_tables/"
count_dz_c = 0
AF_dz_c = 0
min_af_dz_c = 100
max_af_dz_c = 0
count_dz_a = 0
AF_dz_a = 0
min_af_dz_a = 100
max_af_dz_a = 0
count_nondz_c = 0
AF_nondz_c = 0
min_af_nondz_c = 100
max_af_nondz_c = 0
count_nondz_a = 0
AF_nondz_a = 0
min_af_nondz_a = 100
max_af_nondz_a = 0
with open(path + "Final_dz_present_not.txt", encoding = "ISO-8859-1") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[10] != "NA" and line[10] != '':
            if line[2] == "y":
                if line[3] == "causative":
                    count_dz_c +=1
                    AF_dz_c += float(line[10])
                    if float(line[10]) > float(max_af_dz_c):
                        max_af_dz_c = line[10]
                    elif float(line[10]) < float(min_af_dz_c):
                        min_af_dz_c = line[10]
                else:
                    count_dz_a +=1
                    AF_dz_a += float(line[10])
                    if float(line[10]) > float(max_af_dz_a):
                        max_af_dz_a = line[10]
                    elif float(line[10]) < float(min_af_dz_a):
                        min_af_dz_a = line[10]
            else:
                if line[3] == "causative":
                    count_nondz_c +=1
                    AF_nondz_c += float(line[10])
                    if float(line[10]) > float(max_af_nondz_c):
                        max_af_nondz_c = line[10]
                    elif float(line[10]) < float(min_af_nondz_c):
                        min_af_nondz_c = line[10]
                else:
                    count_nondz_a +=1
                    AF_nondz_a += float(line[10]) 
                    if float(line[10]) > float(max_af_nondz_a):
                        max_af_nondz_a = line[10]
                    elif float(line[10]) < float(min_af_nondz_a):
                        min_af_nondz_a = line[10]
print(count_dz_c)
print(AF_dz_c/count_dz_c)
print(min_af_dz_c)
print(max_af_dz_c)
print(count_dz_a)
print(AF_dz_a/count_dz_a)
print(min_af_dz_a)
print(max_af_dz_a)
print(count_nondz_c)
print(AF_nondz_c/count_nondz_c)
print(min_af_nondz_c)
print(max_af_nondz_c)
print(count_nondz_a)
print(AF_nondz_a/count_nondz_a)        
print(min_af_nondz_a)
print(max_af_nondz_a)
