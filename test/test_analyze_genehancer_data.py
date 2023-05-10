import pandas
import numpy as np
import os
import sys
import time

from line_profiler import LineProfiler #type:ignore

sys.path.append('../src/regelem')
import analyze_genehancer as agd
import files


# * Stuff needs to be defined here for LineProfiler. Else it wont find variables in scope. 
input_folder = "../input/genehancer/GeneHancer_v5.14.sorted.gff"
folder = "./output/genehancer_dataframe/GeneHancer_res1000.csv"

test_folder = "./output/genehancer_dataframe/line_profiler_test.csv"

df = agd.read_genehancer_to_df(input_folder)

df = df[df['#chrom'] == "chr1"]

prom_df = df[df['feature name'] == 'Promoter']
enh_df = df[df['feature name'] == 'Enhancer']

def main():
    start_time = time.time()
    dataframe = agd.genehancer_data_to_bins_multiprocess(input_folder,1000,folder)
    end_time = time.time() - start_time

    print(dataframe)
    print(f"genehancer_data_to_bins_multiprocess finished in {end_time}")

    # #test_folder = "./output/genehancer_dataframe/CHRY_line_profiler.csv"

    # df = agd.read_genehancer_to_df(input_folder)

    # df = df[df['#chrom'] == "chr1"]

    # prom_df = df[df['feature name'] == 'Promoter']
    # enh_df = df[df['feature name'] == 'Enhancer']

    # print(prom_df)
    # print(enh_df)
    # #agd._analyze(prom_df,enh_df,5000,test_folder)

    # lp = LineProfiler()
    # lp.add_function(agd._analyze)
    # lp.run('agd._analyze(prom_df,enh_df,5000,test_folder)')
    # lp.print_stats()

    #exit()

    
    # lp = LineProfiler()
    # lp.add_function(agd.genehancer_data_to_bins_multiprocess)
    # lp.run('agd.genehancer_data_to_bins_multiprocess(input_folder,5000,folder,target_chrom="chrY")')
    # lp.print_stats()


    # genehancer_df = agd.read_genehancer_to_df(input_folder)

    # gene_enhancer_hash = agd.make_gene_enhancer_hash(genehancer_df)

    # agd.gather_temp_files(folder)

    pass


if __name__ == "__main__":
    main()