import numpy as np
import pandas as pd
import sys

sys.path.insert(0,"../src/regelem")
import files
import compare_genehancer_to_cooler as cgtc

cooler_input = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.bed"
genehancer_input = "../input/dataframes/genehancer/GeneHancer_res5000.csv"

genehancer_gff_input = "../input/genehancer/GeneHancer_v5.14.gff"
pls_end_input = "../input/GRCh38-cCREs.PLSELS.bed"

def main():

    distance_limit = 5000

    cgtc.stolen_code_from_roberto(pls_end_input,genehancer_gff_input)

    cgtc.assosiate_encode_genehancer_id(genehancer_gff_input,pls_end_input,distance_limit)


    #cgtc.compare(genehancer_input,cooler_input,distance_limit)

    

if __name__ == "__main__":
    main()