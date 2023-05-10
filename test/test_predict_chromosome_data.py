import sys
import os
from os.path import dirname, basename, isfile, join

sys.path.insert(0, '..')
from Constants import *
import predict_chromosome_data
#from predict_chromosome_data import *

def main():
    test_prediction_of_promoter_enhancer()

    print(pixels_on_diagonal(10,3))
    print(pixels_on_diagonal(2_489_564,30_000))

def pixels_on_diagonal(chrom_size : int, diagonal_width: int):
    n_of_pixels = 0

    for i in range(chrom_size-diagonal_width):
            n_of_pixels += 2 * diagonal_width + 1
    for i in range(diagonal_width):
            n_of_pixels += 1 + i * 2
    
    return n_of_pixels

def test_prediction_of_promoter_enhancer():
    start_time = time.time()

    bed_file_path = input_folder + "H1-hESC.7group.bed"
    #bed_file_path = input_folder + "test_bedfile_multiple_bins.bed"
    chrom_sizes_file_path = input_folder + "hg38.chrom.sizes"

    # Chrom name present in test file
    chrom_name = "chr1"

    # Read bed file
    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])
    # Get the size of the chromosome

    chrom_size = int(chrom_sizes_dataframe[chrom_sizes_dataframe['chrom'] == chrom_name]['size'])
    print(chrom_size)

    # Hard code a fake chrom size for testing purposes
    #chrom_size = 1_100_000
    resolution = 50

    promoter_enhancer_dataframe = pd.read_csv(bed_file_path,delim_whitespace=True,header=None,names=DATAFRAME_COLUMNS_BED)
    # Filter 'type' column to only include PLS, pELS or dELS (makes it easier to check this column later)
    promoter_enhancer_dataframe = pelg.filter_type_in_dataframe(promoter_enhancer_dataframe)
    # Split dataframe into one for promoters and one for enhancers
    promoter_dataframe, enhancer_dataframe = pelg.split_df_to_pls_els(promoter_enhancer_dataframe)

    predict_chromosome_data.calc_reg_pairs_multiprocess(   promoter_dataframe=promoter_dataframe,enhancer_dataframe=enhancer_dataframe,
                                        chrom_name=chrom_name,chrom_size=chrom_size, resolution=resolution, 
                                        output_path = f'./output/predicted_chromosome_data_test_{chrom_name}_{resolution}.csv')

    end_time = time.time() - start_time

    print(f"Finished in {end_time} seconds")

if __name__ == "__main__":
    main()