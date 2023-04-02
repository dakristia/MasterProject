
import numpy as np
import pandas as pd
import os
import sys
import Constants

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

def save_dataframe(dataframe : pd.DataFrame, file_path : str, numpy_columns : np.array = False):
    """Save pandas dataframe to file

    Args:
        dataframe (pd.DataFrame): _description_
        filepath (str): _description_
    """
    print("Saving dataframe to filename:", file_path)

    # * Create path and save dataframe to csv
    create_dir_for_path(file_path)

    dataframe = dataframe.copy()

    if numpy_columns:
        index = 0
        for column_name in numpy_columns:
            print(f"Saving {column_name} as numpy array")
            column = dataframe[column_name]

            column_array = np.array([column])[0]

            new_folder = os.path.splitext(file_path)[0]
            numpy_path = f'{new_folder}/numpies/{column_name}.npy'
            create_dir_for_path(numpy_path)
            np.save(numpy_path,column_array)

            # * Replace elements with filepath of array instead. 
            column_file_name_array = np.array([numpy_path for _ in column_array])
            dataframe[column_name] = column_file_name_array

    dataframe.to_csv(file_path)
    
    # * Check if file saved successfully
    check_if_saved = os.path.isfile(file_path)
    if not check_if_saved:
        print(f'Failed to save dataframe to file: {file_path} \n{dataframe}')
    else:
        print(f'Successfully saved dataframe to file at {file_path} (or atleast the file now exists)')

def load_dataframe(filepath : str, numpy_columns : list = False, allow_pickle = False):
    """Attempts to load a dataframe from a file. Also drops "Unnamed: 0" column if present.
    #TODO: Fix so Unnamed: 0 will never occur. This happens because its the old index column.


    Args:
        filepath (): path to dataframe

    Returns:
        pd.Dataframe: The loaded dataframe. False if file isn't found.
    """
    try:
        print("Retrieving dataframe from:", filepath)
        dataframe = pd.read_csv(filepath,index_col=0)
        if "Unnamed:0" in dataframe.columns: dataframe = dataframe.drop("Unnamed: 0", axis=1)
        print("Successfully loaded dataframe from:", filepath)

        index = 0

        if numpy_columns:
            print("Numpy columns parameter True. Attempting to load numpy arrays")
            for column_name in numpy_columns:
                column = dataframe[column_name]
                numpy_path = column[0]
                print(f"Loading: {numpy_path}")
                try:
                    dataframe[column_name] = np.load(numpy_path, allow_pickle=allow_pickle)
                except ValueError as e:
                    e.__str__
                    if str(e) == "Object arrays cannot be loaded when allow_pickle=False":
                        print(f"ENCOUNTERED ERROR: {e}\n Set allow_pickle paramter to true if you trust the source.")
                        return False
                    else:
                        print(e)

        return dataframe

    except FileNotFoundError:
        print("Failed to load dataframe from:", filepath)
        return False

def extract_pls_els_from_bed(bed_file_name: str, split : bool = False) -> pd.core.frame.DataFrame:
    """Read bed file and extract PLS and ELS row data.

    Args:
        bed_file_name (str): path to bed file.
        split (bool, optional): Wether to split the data into two dataframes. One for PLS and one for ELS. If True, returns tuple of dataframes. Defaults to False.

    Returns:
        pd.core.frame.DataFrame: Dataframe with PLS and ELS.
    """
    print("Creating dataframe with promoters and enhancers from .bed file")
    df : pd.core.frame.DataFrame = pd.read_csv(bed_file_name, delim_whitespace=True, header = None, names = Constants.DATAFRAME_COLUMNS_BED)

    pls_dataframe = df.loc[df["type"].str.contains("PLS")]
    els_dataframe = df.loc[df["type"].str.contains("ELS")]
    
    if split: return pls_dataframe, els_dataframe

    df = pd.concat([pls_dataframe,els_dataframe],ignore_index=True).drop_duplicates().sort_values(by=['chrom','chromStart','chromEnd'])
    return df