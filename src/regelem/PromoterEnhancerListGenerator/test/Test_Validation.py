import sys
import pytest  # type:ignore
sys.path.insert(0, '..')
from File_Functions import *
from Functions import *

test_folder = "./testInput/"


def compare_modle_mcool():
    """ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"""

    chrom_name = "chr19"

    bed_file = test_folder + "ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"

    # mcool validation data

    filename_mcool = test_folder + "4DNFI9GMP2J8.mcool::/resolutions/5000"
    mcool_file = read_mcool_file(filename_mcool)

    mcool_2Dselector = mcool_file.matrix(
        balance=False, as_pixels=True, join=True)

    #mcool_matrix = mcool_2Dselector.fetch((chrom_name,100_000,1_000_000))
    
    
    mcool_interaction_dataframe = create_list_of_chrom(filename_mcool, bed_file, chrom_name)
    
    # Modle test data
    filename_modle = test_folder + "test.cool"
    modle_cool_file, modle_2Dselector = read_cool_file(filename_modle)

    modle_2Dselector = modle_cool_file.matrix(
        balance=False, as_pixels=True, join=True)
    
   # modle_matrix = modle_2Dselector.fetch((chrom_name,100_000,1_000_000))

    modle_interaction_dataframe = create_list_of_chrom(filename_modle, bed_file, chrom_name)

    ## Data Comparison

    # print(mcool_matrix)
    # print(modle_matrix)



    print(mcool_interaction_dataframe)
    print(modle_interaction_dataframe)


def main():
    compare_modle_mcool()

main()