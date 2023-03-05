import sys
import os
from os.path import dirname, basename, isfile, join

sys.path.insert(0, '..')
from Constants import *
from predict_chromosome_data import *


def test_prediction_of_promoter_enhancer():
    #bed_file_path = input_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"
    bed_file_path = input_folder + "test_bedfile_multiple_bins.bed"
    chrom_sizes_file_path = input_folder + "hg38.chrom.sizes"

    # Chrom name present in test file
    chrom_name = "chr1"

    # Read bed file
    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])
    # Get the size of the chromosome
    chrom_size = int(chrom_sizes_dataframe[chrom_sizes_dataframe['chrom'] == chrom_name]['size'])
    
    # Hard code a fake chrom size for testing purposes
    chrom_size = 1_100_000
    resolution = 100_000

    promoter_enhancer_dataframe = pd.read_csv(bed_file_path,delim_whitespace=True,header=None,names=DATAFRAME_COLUMNS_BED)
    # Filter 'type' column to only include PLS, pELS or dELS (makes it easier to check this column later)
    promoter_enhancer_dataframe = pelg.filter_type_in_dataframe(promoter_enhancer_dataframe)
    # Split dataframe into one for promoters and one for enhancers
    promoter_dataframe, enhancer_dataframe = pelg.split_df_to_pls_els(promoter_enhancer_dataframe)

    calculate_promoter_enhancer_bins(   promoter_dataframe=promoter_dataframe,enhancer_dataframe=enhancer_dataframe,
                                        chrom_name=chrom_name,chrom_size=chrom_size, resolution=resolution)
