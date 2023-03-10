
import numpy as np
import pandas as pd
import os
import sys
from Constants import *

def create_dir_for_path(path : str):
    """Create folder along path to target if they don't already exist.
    Note: If end of path is a file, it will not create the file itself. Only the folders that lead to it

    Args:
        path (str): path to create folders for.
    """

    path_delim = "/"

    path_array = path.split(path_delim)

    cur_path = ""
    if path[0] == path_delim:
        cur_path+=path_delim

    for i in range(len(path_array)):
        cur_dir = path_array[i]
        cur_path += cur_dir

        # * If not at a relative path (means dir already exists)
        if cur_dir != "." and cur_dir != "..":
            # * If not at the end of path and if dir doesn't exist 
            if not len(cur_path) == len(path) and not os.path.isdir(cur_path):
                    print("Creating directory at:", cur_path)
                    # * Create dir 
                    os.mkdir(cur_path)
        cur_path += path_delim
        print(cur_path)

def save_matrix(matrix : np.ndarray, file_path : str):
    """Save numpy matrix to file.

    Args:
        matrix (np.ndarray): matrix to save
        file_path (str): path to save location
    """
    print("Saving matrix to:", file_path)
    create_dir_for_path(file_path)
    np.save(file_path, matrix)

def load_matrix(file_path: str) -> np.array:
    """Load numpy matrix from file

    Args:
        file_path (str): path to file

    Returns:
        np.array: the loaded matrix
    """
    try:
        print("Retrieving matrix from:", file_path)
        matrix = np.load(file_path,allow_pickle=True) #! TODO: allow_pickle is a security issue
        print("Successfully loaded matrix from:", file_path)
        return matrix 
    except OSError:
        print("Failed to load matrix from:",file_path)
        return False

def save_dataframe(dataframe : pd.DataFrame, file_path : str):
    """Save pandas dataframe to file

    Args:
        dataframe (pd.DataFrame): _description_
        filepath (str): _description_
    """
    print("Saving dataframe to filename:", file_path)

    # * Create path and save dataframe to csv
    create_dir_for_path(file_path)
    dataframe.to_csv(file_path)
    
    # * Check if file saved successfully
    check_if_saved = os.path.isfile(file_path)
    if not check_if_saved:
        print(f'Failed to save dataframe to file: {file_path} \n{dataframe}')

def load_dataframe(filepath : str):
    try:
        print("Retrieving dataframe from:", filepath)
        dataframe = pd.read_csv(filepath)
        print("Successfully loaded dataframe from:", filepath)
        return dataframe
    except FileNotFoundError:
        print("Failed to load dataframe from:", filepath)
        return False

def extract_pls_els_from_bed(bed_file_name: str, split : bool = False) -> pd.core.frame.DataFrame:
    print("Creating dataframe with promoters and enhancers from .bed file")
    df : pd.core.frame.DataFrame = pd.read_csv(bed_file_name, delim_whitespace=True, header = None, names = DATAFRAME_COLUMNS_BED)

    pls_dataframe = df.loc[df["type"].str.contains("PLS")]
    els_dataframe = df.loc[df["type"].str.contains("ELS")]
    
    if split: return pls_dataframe, els_dataframe

    df = pd.concat([pls_dataframe,els_dataframe],ignore_index=True).drop_duplicates().sort_values(by=['chrom','chromStart','chromEnd'])
    return df