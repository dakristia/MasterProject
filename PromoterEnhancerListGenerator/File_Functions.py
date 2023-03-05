from Constants import *
from numba import njit
import re  # Regular Expressions
import cooler  # type:ignore
import pandas as pd  # type:ignore
import pandas.core.frame as pdFrame
import numpy as np
import subprocess



def grep_bed(filename: str):
    bashcommand = ["echo","'PLS\|ELS'", filename]
    print(bashcommand)
    #bashcommand = "grep 'PLS\|ELS' " + filename#, ">", filename.replace(".bed","_grepped.bed")]
    print(subprocess.run(bashcommand, capture_output=True).stdout)

# Reads a .bed file and adds its content to a 2D-list
def read_Bedfile(filename: str) -> list:
    """Read the contents of a .bed file, save each line as element in a list, return list

    Keyword arguments:
    filename -- a filename string to a .bed file
    """

    bedContentList: list = []
    with open(filename) as file:
        for line in file:
            bedContentList.append(line.strip().split())

    return bedContentList


def read_Gtfile(filename: str) -> list:
    """Read the contents of a .gtf file 

    Keyword arguemnts:
    filename -- a filename string to a .gtf file
    """

    gtfContent: list = []
    with open(filename) as file:
        for line in file:
            gtfContent.append(line.strip().split())

    return gtfContent

def read_mcool_file(mcool_file_path: str):
    modle_cool: cooler.api.Cooler = cooler.Cooler(mcool_file_path)
    print(modle_cool)
    print(modle_cool.chromnames)

    return modle_cool
    #print(modle_cool.matrix().fetch("chr19"))
    #TODO first manually compare this data to the stuff produced by modle
    #TODO we need to focus on a particular region. Any regions which aren't noted by modle are irrelevant

def read_cool_file(cooler_file_path: str) -> tuple:
    """Reads the info in a .cool file and prints it to the terminal. 

    Keyword arguments:
    cooler_file_path -- a filepath to a .cool file

    Returns:
    A tuple containing the cooler object of the file, and a 2D selector for the matrix.
    """
    # Attempt to import data with cooler
    modle_cool: cooler.api.Cooler = cooler.Cooler(cooler_file_path)
    resolution: np.int64 = modle_cool.binsize
    chromnames: list = modle_cool.chromnames
    chromsizes: pd.core.series.Series = modle_cool.chromsizes
    info: dict = modle_cool.info

    print("Cooler stats -----------------------")
    print("Resolution:", resolution, "\n\n")
    print("chromnames:", chromnames, "\n\n")
    print("chromssizes:", chromsizes, "\n\n")
    print("info:", info, "\n\n")

    # This is a 2D selector object, very light weight
    cooler_2Dselector = modle_cool.matrix(
        balance=False, as_pixels=True, join=True)

    return (modle_cool, cooler_2Dselector)


def extract_pls_els_from_bed(bed_file_name: str) -> pd.core.frame.DataFrame:
    print("Creating dataframe with promoters and enhancers from .bed file")
    df : pd.core.frame.DataFrame = pd.read_csv(bed_file_name, delim_whitespace=True, header = None, names = DATAFRAME_COLUMNS_BED)
    return df


def extract_Extrusion_Barrier_From_Bed(bedList: list) -> dict:
    """Extract a dict of all lines from a bedFile


    bedList -- a list containing strings, where each element is a line in a .bed file
            -- a filename string to a .bed file
    """

    if type(bedList) is str:
        bedList = read_Bedfile(bedList)

    extrusion_dict = {}

    for bedLine in bedList:
        chromName = bedLine[BED_CHROM]
        if chromName not in extrusion_dict.keys():
            extrusion_dict[chromName] = []

        extrusion_dict[chromName].append(bedLine)

    return extrusion_dict


def write_to_file_bedpe(dataframe: pd.DataFrame, out_file_path: str) -> None:
    """
    Writes the contents of a dataframe to a .bedpe file. Note that this doesn't include validation data for now

    Format:
    chrom1 prom_start prom_end chrom2 enh_start enh_end name, modle_count, strand1, strand2

    Arguments:
        dataframe: a pandas dataframe. Requires the columns: 'chrom', 'p_name', 'e_name' , 'count', 'freq'
        out_file_path: the file path of the output file
    Returns:
        None

    """

    contact_list = []

    chrom_index = DATAFRAME_COLUMNS_INTERNAL.index("chrom") + 1
    enh_start_index = DATAFRAME_COLUMNS_INTERNAL.index("enh_start") + 1
    enh_end_index = DATAFRAME_COLUMNS_INTERNAL.index("enh_end") + 1
    prom_start_index = DATAFRAME_COLUMNS_INTERNAL.index("prom_start") + 1
    prom_end_index = DATAFRAME_COLUMNS_INTERNAL.index("prom_end") + 1
    enh_start_index = DATAFRAME_COLUMNS_INTERNAL.index("enh_start") + 1
    enh_start_index = DATAFRAME_COLUMNS_INTERNAL.index("enh_end") + 1
    prom_name_index = DATAFRAME_COLUMNS_INTERNAL.index("prom_name") + 1
    enh_name_index = DATAFRAME_COLUMNS_INTERNAL.index("enh_name") + 1 
    count_index = DATAFRAME_COLUMNS_INTERNAL.index("modle_count") + 1

    for dataframe_row in dataframe.itertuples():
        
        contact_list.append([dataframe_row[chrom_index],
                            dataframe_row[prom_start_index], dataframe_row[prom_end_index],
                            dataframe_row[chrom_index],
                            dataframe_row[enh_start_index],dataframe_row[enh_end_index],
                            str(dataframe_row[prom_name_index]) + "(prom)-(enh)" + str(dataframe_row[enh_name_index]),
                            dataframe_row[count_index],
                            ".",
                            "."])
    
    
    with open(out_file_path, 'w') as out_file:
        for line in contact_list:
            for element in line:
                out_file.write(str(element) + " ")
            out_file.write("\n")
    
    return