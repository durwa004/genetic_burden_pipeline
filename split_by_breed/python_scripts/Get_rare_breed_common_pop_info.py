import gzip
import os

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect"
for filename in os.listdir(directory):
    if filename.endswith("_rare_population.vcf.gz"):
        with gzip.open(directory + filename, "rt") as input_file, open("breed_common_rare_population_summary.txt", "w") as output_file:
            a = filename.split("_")
            breed = a[1]
            for line in input_file:
                line = line.rstrip)"\n").split("\t")
                if "#" in line[0]:
                    next
                else:
                    


        if filename == "thesis_intersect.vcf.gz":
            next
        
        else:
            with gzip.open("thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file, open("SnpEff_variant_type_all.txt", "w") as output_file:
    print("#CHROM\tPOS\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            next
        else:
            if "NC_" in line[0]:
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                AC = bc[1]
                cd = ab[1].split("AF=")
                AF = cd[1]
                de = line[7].split("ANN=")
                bc = de[1].split("|")
                consequence = bc[1]
                impact = bc[2]
                gene = bc[4]
                if "splice" in consequence:
                    consequence = "splice_region_variant"
                elif "synonymous" in consequence or "stop" in consequence or "start" in consequence or "missense" in consequence or "initiator" in consequence or "frameshift" in consequence or "inframe" in consequence:
                    consequence = "exonic"
                elif "gene_fusion" in consequence:
                    consequence = "gene_fusion"
                print(line[0], line[1], AC, AF, consequence, impact, gene, file = output_file, sep = "\t")
            else:
                next