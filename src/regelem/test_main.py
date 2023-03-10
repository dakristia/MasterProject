import files

import predict_chromosome_data
import calculate_general_statistics_of_cool_data
import PromoterEnhancerListGenerator as pelg


import os
import pandas as pd
import time

def main():
    #save_matrix()
    #create_dir_for_path("./output/lmao/haha")
    chrom_name = "chr19"

    cooler_file_path = "../../input/" + "outfile_binsize5000_tcd7.31.cool"
    bed_file_path = "../../input/" + "H1-hESC.7group.bed"
    chrom_sizes_file_path = "../../input/" + "hg38.chrom.sizes"
    
    if os.path.isfile("../output/test_dataframe_DELETEME.csv"):
        interaction_dataframe = files.load_dataframe("../output/test_dataframe_DELETEME.csv")
    else:
        interaction_dataframe = pelg.Dataframe_Functions.note_promoter_enhancer_interactions_multiprocess(cooler_file_path,bed_file_path,chrom_name)
        files.save_dataframe(interaction_dataframe, "../output/test_dataframe_DELETEME.csv")


    calculate_general_statistics_of_cool_data.calculate_average_counte_per_regelem_bin(interaction_dataframe)
    exit()


    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])
    chrom_size = int(chrom_sizes_dataframe[chrom_sizes_dataframe['chrom'] == chrom_name]['size'])
    resolution = 5000

    start_time = time.time()

    predict_chromosome_data.calculate_promoter_enhancer_bins_multiprocess("../../input/" + "H1-hESC.7group.bed", chrom_name, chrom_size, 5000, 10, 20, f"./output/predicted_regulatoryinteraction_{chrom_name}_{resolution}")

    end_time = time.time() - start_time
    print("finished in", end_time,"seconds")
    pass

if __name__ == "__main__":
    main()