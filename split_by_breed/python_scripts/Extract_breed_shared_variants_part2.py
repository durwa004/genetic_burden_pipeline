import argparse
import os
import gzip

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir with breed thesis intersect files  [required]")
    parser.add_argument(
            "-i", "--inp",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Filename of output of Extract_breed_shared_variants.py   [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input1 = args.inp

    #Want breed and number of variants shared by each breed
    breeds = {}
    for filename in os.listdir(data):
         if filename.endswith(".list"):
             a = filename.split("_")
             if len(a) <3:
                 breeds[a[0]] = {}

    for breed in breeds.keys():
        with open(input1, "r") as input_file:
            for line in input_file:
                line = line.rstrip("\n").split("\t")
                line1 = line[2:]
                for i in range(len(line1)):
                    if line1[i] == breed:
                        for i in range(len(line1)):
                            if line1[i] != breed:
                                if line1[i] in breeds[breed].keys():
                                    a = int(breeds[breed][line1[i]]) + 1
                                    breeds[breed][line1[i]] = a
                                else:
                                    breeds[breed][line1[i]] = 1

    breed_l = list(set(breeds.keys()))
    with open(data + "/breed_breed_shared_unique_variants.txt", "w") as output_file:
        print("breed", "\t".join(breed_l), sep = "\t", file = output_file)
        for key in breeds.keys():
            n_variants =  list(set(breeds.keys()))
            print(key)
            for key1 in breeds[key].keys():
                for i in range(len(breed_l)):
                    if key1 == breed_l[i]:
                        n_variants[i] = breeds[key][key1]
                for i,v in enumerate(n_variants):
                    if isinstance(v,int) == True:
                        next
                    else:
                        n_variants[i] = "NA"
            ab = [str(i) for i in n_variants]
            print(key, "\t".join(ab), sep = "\t", file = output_file) 
