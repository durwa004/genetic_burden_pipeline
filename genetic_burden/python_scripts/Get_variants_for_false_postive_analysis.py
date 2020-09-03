import random

gb = {}
lof = {}
with open("/Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/genetic_burden_details_brief_tidy.txt", "r") as input_file:
    header = input_file.readline()
    header = header.rstrip("\n").split("\t")
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2]
        if a in gb.keys():
            print(a)
        else:
            gb[a] = line[6]

with open("/Users/durwa004/Documents/Postdoc/Projects/LOF_variants/lof_tidy_brief.txt", "r") as input_file:
    header = input_file.readline()
    header = header.rstrip("\n").split("\t")
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[4] + ":" + line[5]
        if a in lof.keys():
            print(a)
        else:
            lof[a] = line[12]

#Get intersect and then ones that are only in GB set or only in LOF set.
intersect = list(set(list(gb.keys())) & set(list(lof.keys())))

gb_only = list(set(list(gb.keys())) - set(intersect))
lof_only = list(set(list(lof.keys())) - set(intersect))
#Get randomly selected 5% of LOF and GB only variants.

len(gb_only) * 0.05
len(intersect) * 0.05
gb_only = list(set(gb_only))
gb_r = random.sample(gb_only, 42)
lof_r = random.sample(intersect, 249)
gb_r = gb_r + lof_r

with open("genetic_burden_details_brief_tidy.txt", "r") as input_file, open("gb_validation_random_set.txt", "w") as f:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i,v in enumerate(gb_r):
                a = v.split(":")
                if line[1] == a[0] and line[2] == a[1]:
                    print("\t".join(line), file = f)
        for i1,v1 in enumerate(lof_only):
            for item in range(len(gb_r)):
                if v1 == gb_r[item]:
                    print("\t".join(line), file = f)

#Get overlap between gb_r and the variants that kelsey has been working on.
kelsey = {}
with open("kelsey_done.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        kelsey[a] = "NA"

to_do = list(set(gb_r) - set(list(kelsey.keys())))

with open("genetic_burden_details_brief_tidy.txt", "r") as input_file, open("gb_validation_to_do.txt", "w") as f:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i,v in enumerate(to_do):
                a = v.split(":")
                if line[1] == a[0] and line[2] == a[1]:
                    print("\t".join(line), file = f)

