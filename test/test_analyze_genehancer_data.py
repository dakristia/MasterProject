import pandas
import numpy as np
import os
import sys
import time

sys.path.append('../src/regelem')
import analyze_genehancer_data as agd
import files

input_folder = "../input/genehancer/GeneHancer_v5.14.sorted.gff"
folder = "./output/genehancer_dataframe/CHR1_again.csv"

def main():
    
    start_time = time.time()
    dataframe = agd.genehancer_data_to_bins_multiprocess(input_folder,5000,folder,target_chrom="chr1")
    end_time = time.time() - start_time
    print("genehancer_data_to_bins_multiprocess finished in end_time")


    genehancer_df = agd.read_genehancer_to_df(input_folder)

    #gene_enhancer_hash = agd.make_gene_enhancer_hash(genehancer_df)

    exit()
    agd.gather_temp_files(folder)

    pass



if __name__ == "__main__":
    main()