#Extract random validation list from variants that have been analyzed by Kelsey and the other diseases

validation = {}
details = {}
with open("gb_validation_random_set.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2]
        validation[a] = "NA"
        details[a] = "NA"

with open("analyzed_variants_brief.txt", "r", encoding="ISO-8859-1") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        if a in validation.keys():
            if line[2] == "Y" or line[2] == "y":
                validation[a] = "y"
            elif line[2] == "N" or line[2] == "n" or line[2] == "n ":
                validation[a] = "n"
                if line[3] == 'no annotated transcriptions' or line[3] == "Doesn't blast to correct human reference gene" or line[3] == 'no significant similarity found' or line[3] == 'no annotated transcribed products' or line[3] == 'shortened and not matching correct gene' or line[3] == 'no significant similarities found' or line[3] == 'shortened and not matching gene' or line[3] == 'best aligned BLAST sequence only about 40bp long':
                    details[a] = "incorrect_annotation"
                elif line[3] == 'exon cut in half' or line[3] == 'last 10bp cut off' or line[3] == 'shortened' or line[3] == 'shortenedÊ' or line[3] == 'first 10 bp cut off' or line[3] == "." or line[3] == 'first part of exon cut off' or line[3] == 'last part of exon cut off' or line[3] == 'had to stitch exons together; first is cut off':
                    details[a] = "incomplete_annotation"
                else:
                   print(line)
            elif line[2] == "?":
                validation[a] = "?"

v = list(validation.values())
print(v.count("n"))
print(v.count("y"))
print(v.count("?"))
print(v.count("NA"))

d = list(details.values())
print(d.count("incorrect_annotation"))
print(d.count("incomplete_annotation"))
