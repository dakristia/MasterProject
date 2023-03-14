import pytest  # type:ignore
import sys
sys.path.insert(0, '..')
from Functions import *
from Dataframe_Functions import *

test_folder = "./testInput/"
test_bed_file = "2022-08-30_cisRegElements_ENCFF428XFI_ENCFF280PUF_ENCFF469WVA_ENCFF644EEX_grepped.7group.bed"
test_cool_file = "OUTPUT.cool"

def test_get_extend_of_PLSELS():
    pass

#! Unfinished
def test_find_pls_els_regions():

    chrom_name = "chr1"

    bed_file = test_folder + test_bed_file

    pls_els_bed_list = read_Bedfile(bed_file)

    #pls_els_dict = extract_PLS_ELS_From_Bed_List()

#! Unfinished
def test_filter_type_in_dataframe():

    bed_file = test_folder + test_bed_file

    cooler_file = test_folder + test_cool_file

    modle_cool, cooler_2Dselector = read_cool_file(cooler_file)
    chrom_names = modle_cool.chromnames

    bed_df = extract_pls_els_from_bed(bed_file)

    dataframe = filter_type_in_dataframe(bed_df)

    #assert dataframe.shape == (5,11)

    #assert list(dataframe.columns) == (DATAFRAME_COLUMNS_BED)

    #assert dataframe['type'].tolist() == ['PLS','PLS','pELS','dELS','PLS']
    
def test_note_interactions():
    bed_file = test_folder + test_bed_file
    cooler_file = test_folder + test_cool_file

    modle_cool, cooler_2Dselector = read_cool_file(cooler_file)

    chrom_names = modle_cool.chromnames
    bed_df = extract_pls_els_from_bed(bed_file)
    dataframe = filter_type_in_dataframe(bed_df)

    chrom_name = chrom_names[0]
    note_interactions(cooler_2Dselector, dataframe, chrom_name)


def test_create_list_of_chrom():
    bed_file = test_folder + test_bed_file
    modle_file_path = test_folder + test_cool_file
    mcool_file_path = test_folder + "4DNFI9GMP2J8.mcool"

    chrom_name = "chr19"
    start = "400_000"
    end = "1_000_000"


    dataframe = create_list_of_chrom(modle_file_path, bed_file, chrom_name, start, end, mcool_file_path)
    print(dataframe)


