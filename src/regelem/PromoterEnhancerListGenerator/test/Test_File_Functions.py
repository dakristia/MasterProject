import sys
import pytest  # type:ignore
sys.path.insert(0, '..')
from File_Functions import *
from Functions import *


test_folder = "./testInput/"


def test_read_Bed():

    bed_file = test_folder + "PLS_pdELS.7_testgroup.bed"

    filename = bed_file
    bed_list = read_Bedfile(filename)

    # Checking if function returns correct type
    assert isinstance(bed_list, list)

    # Checking if list has correct shape and length
    assert len(bed_list) == 6
    assert len(bed_list[0]) == 11

    # Hard checking the first line in the file
    assert bed_list[0][0] == 'chr1'
    assert bed_list[0][1] == '725000'
    assert bed_list[0][2] == '725300'
    assert bed_list[0][3] == 'EH38E2000000'
    assert bed_list[0][4] == '0'
    assert bed_list[0][5] == '.'
    assert bed_list[0][6] == '725000'
    assert bed_list[0][7] == '725346'
    assert bed_list[0][8] == '255,0,0'
    assert bed_list[0][9] == 'PLS,CTCF-bound'
    assert bed_list[0][10] == 'All-data/Full-classification'

    # Makes sure we don't have more entries than we should
    with pytest.raises(IndexError):
        assert bed_list[0][11]

def test_grep_bed():
    bed_file = test_folder + "ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    bed_file = "temptext.txt"
    grep_bed(bed_file)


def test_extract_pls_els_from_bed():
    bed_file = test_folder + "PLS_pdELS.7_testgroup.bed"
    extract_pls_els_from_bed(bed_file)


def test_read_Gtf():
    filename = test_folder + "hg38.refGene.gtf"
    gtf_list = read_Gtfile(filename)

    assert type(gtf_list) == list
    assert len(gtf_list) > 0
    # This will probably need to change as we update how we read gtf
    assert len(gtf_list[0]) == 14

    assert True


def test_read_cool_file():
    filename = test_folder + "test.cool"

    try:
        cooler_object, selector = read_cool_file(filename)
    except ValueError as exc:
        assert False, f"'Calling and returning read_Modle raises exception:  {exc}"

    assert isinstance(cooler_object, cooler.api.Cooler)
    assert len(cooler_object.chromnames) > 0
    assert len(cooler_object.chromnames) == 22

    assert isinstance(selector, cooler.core.RangeSelector2D)
    assert len(selector) == 578016

    chrom_name = cooler_object.chromnames[0]
    dataframe = selector.fetch((chrom_name, 0, 1_000_000))

    assert dataframe['chrom1'][0] in cooler_object.chromnames



def test_read_mcool_file():
    filename = test_folder + "4DNFI9GMP2J8.mcool::/resolutions/5000"
    cooler_file = read_mcool_file(filename)
    filename = test_folder + "4DNFI9GMP2J8.mcool::/resolutions/2000"
    cooler_file = read_mcool_file(filename)
    # try:
    # except ValueError as exc:
    #     assert False, f"'Calling and returning read_Modle raises exception:  {exc}"


def test_mcool_data():
    filename = test_folder + "4DNFI9GMP2J8.mcool::/resolutions/5000"
    cooler_file = read_mcool_file(filename)

    cooler_2Dselector = cooler_file.matrix(
        balance=False, as_pixels=True, join=True)

    data = cooler_2Dselector.fetch(("chr19",100_000,1_000_000))
    print(data)



def test_write_to_file_bedpe():
    bed_file = test_folder + "PLS_pdELS.7_testgroup_big.bed"
    cooler_file = test_folder + "test.cool"
    modle_cool, cooler_2Dselector = read_cool_file(cooler_file)
    chrom_names = modle_cool.chromnames
    bed_df = extract_pls_els_from_bed(bed_file)
    dataframe = filter_type_in_dataframe(bed_df)
    interaction_df = note_interactions(cooler_2Dselector, dataframe, chrom_names)

    write_to_file_bedpe(interaction_df, "test_output.bedpe")

def main():
    pass
    #test_write_to_file_bedpe()
    #test_read_mcool_file()
    #test_mcool_data()

main()