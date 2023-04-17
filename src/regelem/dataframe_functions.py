import os
import warnings
import cooler                   #type:ignore
import matplotlib as mpl        #type:ignore
import matplotlib.pyplot as plt #type:ignore
import numpy as np              #type:ignore
import pandas as pd             #type:ignore
import sys 
import concurrent.futures

import Constants



# Helper method
def split_df_to_pls_els(promoter_enhancer_dataframe : pd.core.frame.DataFrame, column_name = "type") -> tuple:
    """Splits dataframe into two new dataframes. One containing promoter-like elements and one containing enhancer-like.
    Requires the 'type' column. 

    Args:
        promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe to split. Requires 'type' column

    Returns:
        tuple: tuple with promoter and enhancer dataframe
    """
    print("Splitting dataframe into promoter and enhancer dataframe. ")

    pls_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[column_name].str.contains("PLS")]
    els_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[column_name].str.contains("ELS")]
    return (pls_dataframe, els_dataframe)



def downscale_dataframe_to_resolution(dataframe: pd.core.frame.DataFrame,
                                        target_resolution : int = 5000) -> pd.DataFrame:
    """ Downscales a dataframe to the target resolution. 

    Args:
        dataframe (pd.core.frame.DataFrame): Dataframe to scale
        target_resolution (int, optional): Resolution to scale to. Defaults to 5000.

    Returns:
        pd.Dataframe: Dataframe scaled to the target resolution
    """

    ## ! Tacky, bad solution
    # * When dataframes are handled inproperly, an Unnamed: 0 column might be added (probably a stray index column).
    # * If this happens, we delete it.
    if 'Unnamed: 0' in dataframe.columns:
        dataframe = dataframe.drop(columns='Unnamed: 0')

    column_set = False


    try:
        if (dataframe.columns == Constants.DATAFRAME_COLUMNS_COOLER).all() and not column_set:
            columns = Constants.DATAFRAME_COLUMNS_COOLER
            bin1_start_label = "start1"
            bin1_end_label = "end1"
            bin2_start_label = "start2"
            bin2_end_label = "end2"
            count_label = "count"
            column_set = True
    except ValueError: # If column shape mismatch
        pass
    try:
        if (dataframe.columns == Constants.DATAFRAME_COLUMNS_INTERNAL).all() and not column_set:
            columns = Constants.DATAFRAME_COLUMNS_INTERNAL
            bin1_start_label = "bin1_start"
            bin1_end_label = "bin1_end"
            bin2_start_label = "bin2_start"
            bin2_end_label = "bin2_end"
            count_label = "modle_count"
            column_set = True
    except ValueError:
        pass
    if not column_set: 
        print("ERROR Function not implemented for these dataframe columns")
        return False

    resolution = abs(dataframe[bin1_start_label][0] - dataframe[bin1_end_label][0])

    if target_resolution % resolution != 0:
        print("ERROR target_resolution not divisible with original resolution.")
        return False
    elif target_resolution == resolution:
        print(f'ERROR dataframe seemingly already has target resolution {target_resolution}.')
        return dataframe

    print(f'Downscaling dataframe from {resolution} resolution to {target_resolution} ...')

    dataframe_with_target_res = pd.core.frame.DataFrame(columns = columns)

        


    ## * Collect all counts within 1000 res dataframe to a 5000 one
    for row in dataframe.itertuples():
        row_index = row[0]
        bin1_start_index = columns.index(bin1_start_label)
        bin1_start = row[bin1_start_index + 1]
        bin1_end_index = columns.index(bin1_end_label)
        bin1_end = row[bin1_end_index + 1]
        
        bin2_start_index = columns.index(bin2_start_label)
        bin2_start = row[bin2_start_index+ 1]
        bin2_end_index = columns.index(bin2_end_label)
        bin2_end = row[bin2_end_index + 1]

        bin1_start_rounded, bin1_end_rounded = round_up_and_down(bin1_start,target_resolution)
        bin2_start_rounded, bin2_end_rounded = round_up_and_down(bin2_start,target_resolution)


        count_index = columns.index(count_label)
        count = row[count_index + 1]

        target_row = dataframe_with_target_res.loc[(dataframe_with_target_res[bin1_start_label]==bin1_start_rounded) & (dataframe_with_target_res[bin2_start_label]==bin2_start_rounded)]

        if target_row.empty:
                        row_list = list(row)
                        row_list[bin1_start_index + 1] = bin1_start_rounded
                        row_list[bin1_end_index + 1] = bin1_end_rounded
                        row_list[bin2_start_index + 1] = bin2_start_rounded
                        row_list[bin2_end_index + 1] = bin2_end_rounded
                        
                        dataframe_with_target_res = pd.concat([dataframe_with_target_res, pd.core.frame.DataFrame([row_list[1:]], columns=columns)]).reset_index(drop=True)
        else:
            indexes = dataframe_with_target_res.index[(dataframe_with_target_res[bin1_start_label]==bin1_start_rounded) & (dataframe_with_target_res[bin2_start_label]==bin2_start_rounded)]
            index = indexes[0]
            dataframe_with_target_res.iloc[index][count_label] += count
    try:
        if (dataframe.columns == DATAFRAME_COLUMNS_COOLER).all(): dataframe_with_target_res = dataframe_with_target_res.drop_duplicates().sort_values(by=['chrom1','start1','end1','chrom2','start2','end2']).reset_index(drop=True)
    except ValueError:
        pass
    try:
        if (dataframe.columns == DATAFRAME_COLUMNS_INTERNAL).all(): dataframe_with_target_res = dataframe_with_target_res.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)
    except ValueError:
        pass

    return dataframe_with_target_res


def filter_type_in_dataframe(pe_df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    """Removes all unnecessary keywords in 'type' column of dataframe (mainly everything that isn't PLS, dELS, pELS)

    Args:
        pe_df (pd.core.frame.DataFrame): dataframe to be filtered

    Returns:
        pd.core.frame.Dataframe: dataframe where the 'type' column only contains whitelisted words
    """

    print("Filtering type column in dataframe. ")

    whitelist = ['PLS', 'dELS', 'pELS']

    # Generate regex from whitelist
    regex = ''
    for word in whitelist:
        regex += word
        regex += '|'
    regex = regex[:-1]

    # Getting a dataframe with true on all rows containing a word in whitelist
    type_bool_df = pe_df[Constants.DATAFRAME_COLUMNS_BED[9]].str.contains(regex)
    
    # Make a dataframe containing only rows containing 'PLS'
    promoter_enhancer_dataframe = pe_df.loc[type_bool_df]

    # We only care about these words in 'type' column

    # Defining some lambdas to get keywords out of 
    find_in_whitelist = lambda param: param in whitelist
    contains_whitelisted_word = lambda param: next(i for i in param if find_in_whitelist(i))
    split_then_contains_whitelisted_word = lambda param: contains_whitelisted_word(param.split(','))

    # Remove any data from 'type' column that is not in the whitelist
    data = promoter_enhancer_dataframe['type'].map(split_then_contains_whitelisted_word)

    data = data.to_frame()

    #? Can we make it so this functions without raising an error?
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        promoter_enhancer_dataframe['type'] = data['type']

    return promoter_enhancer_dataframe


def check_if_column_sorted (dataframe : pd.DataFrame,
                            col_name : str):
        column_as_numpy = np.array(dataframe[col_name])

        is_sorted = lambda a: np.all(a[:-1] <= a[1:])

        return is_sorted(column_as_numpy)
