import gzip

#Note - will need to do this for both SNPs and Indels
#Goal is to convert the VEP output to a useable text file for analysis.
#Remove chrUn variants too
with gzip.open("/panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype.goldenPath.vep.vcf.gz", "rt") as input_file, open("VEP_variant_type_all.txt", "w") as output_file:
    print("#CHROM\tPOS\tREF\tALT\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            next
        else:
            details = {}
            chr_pos = line[0] + ":" + line[1] + ":" + line[3]
            details[chr_pos] = {}
            if "," in line[4]:
                al = line[4].split(",")
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                bc = bc[1].split(",")
                cd = ab[1].split("AF=")
                cd = cd[1].split(",")
                for i in range(len(al)): #ref:AC:AF
                    xv = al[i] + ":" + bc[i] + ":" + cd[i]
                    details[chr_pos][xv] = {}
                    details[chr_pos][xv][al[i]] = {} 
                    ef = line[7].split("CSQ=")
                    ef = ef[1].split(",")
                    gh = {}
                    for i in range(len(ef)):
                        fg = ef[i].split("|")
                        if ef[i] in gh.keys():
                             hi = fg[1] + ":" + fg[2] + ":"+ fg[3] + ":" + fg[5] + ":" + fg[9]
                             jk = gh[ef[i]] + "," + hi
                             gh[ef[i]] = jk
                        else:
                             hi = fg[1] + ":" + fg[2] + ":" + fg[3] + ":" + fg[5] + ":" + fg[9]
                             gh[ef[i]] = hi
            else:
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                cd = ab[1].split("AF=")
                ef = line[7].split("CSQ=")
                ef = ef[1].split(",")
                gh = {}
                for i in range(len(ef)):
                    fg = ef[i].split("|")
                    if fg[0] in gh.keys():
                         hi = fg[1] + ":" + fg[2] + ":"+ fg[3] + ":" + fg[5] + ":" + fg[9]
                         jk = gh[ef[i]] + "," + hi
                         gh[ef[i]] = jk
                    else:
                         hi = fg[1] + ":" + fg[2] + ":" + fg[3] + ":" + fg[5] + ":" + fg[9]
                         gh[ef[i]] = hi
                for key in gh.keys():
                    kl = gh[key].split(",")
                    for i in range(len(kl)):
                        lm = kl[i].split(":")
                        print(line[0], "\t", line[1], "\t", line[3], "\t", line[4], "\t", bc[1], "\t", cd[1], "\t".join(lm), file = output_f, sep = "")
                       
##Start from here - add in each of the consequences based on whatever the ref is
                de = line[7].split("CSQ=")
                de = de.split(",")
                for i in range(len(de)):
                    bc = de[1].split("|")
                    ALT.append(bc[0])
                    consequence.append(bc[1])
                    impact.append(bc[2])
                for i in range(len(ALT)):
                    print(line[0], line[1], line[3], ALT[i], AC[i], AF[i], consequence, impact, gene,lof, file = output_file, sep = "\t")
            else:
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                AC = bc[1]
                cd = ab[1].split("AF=")
                AF = cd[1]
                de = line[7].split("CSQ=")
                bc = de[1].split("|")
                consequence = bc[1]
                impact = bc[2]
                gene = bc[4]
                if consequence == "splice_acceptor_variant" or consequence == "splice_donor_variant" or consequence == "splice_region_variant" or consequence == "splice_donor_5th_base_variant" or consequence == "splice_donor_region_variant" or "splice_polypyrimidine_tract_variant":
                    consequence = "splice_region"
                elif consequence == "inframe_insertion" or consequence == "inframe_deletion":
                    consequence = "inframe"
                elif "synonymous_variant" == consequence:
                    consequence = "synonymous"
                if "missense_variant" == consequence:
                    consequence = "missense"
                elif "5_prime_UTR_variant" == consequence or consequence == "3_prime_UTR_variant":
                    consequence = "5/3_prime_UTR"
                elif "frameshift_variant" == consequence:
                    consequence = "frameshift"
                elif "start_lost" == consequence:
                    consequence = "start_lost"
                elif "stop_gained" == consequence:
                    consequence = "stop_gained"
                elif "stop_lost" == consequence:
                    consequence = "stop_lost"
                elif "start_retained_variant" == consequence:
                    consequence = "start_retained"
                elif "stop_retained_variant" == consequence:
                    consequence = "stop_retained"
                elif "gene_fusion" in consequence:
                    consequence = "gene_fusion"
                elif "intergenic_variant" == consequence:
                     consequence = "intergenic"
                elif "intron_variant" == consequence:
                     consequence = "intronic"
                elif "upstream_gene_variant" == consequence or "downstream_gene_variant" == consequence:
                     consequence = "Up/downstream"
                else:
                     print(consequence)
                for i in range(len(ab)):
                    ef = ab[i].split("LOF=")
                    if len(ef) > 1:
                        lof = ab[i]
                    else:
                        lof = "n"
                print(line[0], line[1], line[3], line[4], AC, AF, consequence, impact, gene,lof, file = output_file, sep = "\t")
            else:
                next
