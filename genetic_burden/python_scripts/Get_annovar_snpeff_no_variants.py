moderate_count = 0
high_count = 0
low_count = 0
with open("annovar_coding_tidy.txt", "r") as input_file, open("annovar_GB_info.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if len(line) > 1:
            if line[10] == "HIGH":
                high_count +=1
            elif line[10] == "MODERATE":
                moderate_count +=1
            elif line[10] == "LOW":
                low_count += 1
    print("HIGH:", high_count, "\nMODERATE:", moderate_count, "\nLOW:", low_count, file = output_file, sep="")
     
moderate_count = 0
high_count = 0
low_count = 0
with open("SnpEff_coding_tidy.txt", "r") as input_file, open("SnpEff_GB_info.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if len(line) > 1:
            if line[10] == "HIGH":
                high_count +=1
            elif line[10] == "MODERATE":
                moderate_count +=1
            elif line[10] == "LOW":
                low_count += 1
    print("HIGH:", high_count, "\nMODERATE:", moderate_count, "\nLOW:", low_count, file = output_file, sep="")
