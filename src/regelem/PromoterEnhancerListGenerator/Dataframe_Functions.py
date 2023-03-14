import os
import threading
import warnings
import cooler                   #type:ignore
import matplotlib as mpl        #type:ignore
import matplotlib.pyplot as plt #type:ignore
import numpy as np              #type:ignore
import pandas as pd             #type:ignore
import sys 
import subprocess
import concurrent.futures

#from os.path import dirname, basename, isfile, join

#from .File_Functions import *
from . import (File_Functions, Functions, General_Functions)
import Constants
from .File_Functions import *
from .Functions import *
from .General_Functions import *
import files

def extract_original_dataframe_from_cool(
                        cool_file_path : str,
                        chrom_name : str,
                        start : int = None,
                        end : int = None
                        ) -> pd.core.frame.DataFrame:
    cooler_object : cooler.api.Cooler = cooler.Cooler(cool_file_path)

    cooler_2D_selector = cooler_object.matrix(balance=False, as_pixels=True, join=True)

    cooler_dataframe = cooler_2D_selector.fetch((chrom_name,start,end))

    return cooler_dataframe

#TODO: Check how many times this is used
def filter_cooler_to_regelement_dataframe(
                        cool_file_path : str,
                        bed_file_path : str,
                        chrom_name : str, 
                        resolution : int = False,
                        start : int = None,
                        end : int = None,
                        max_memory : int = None) -> pd.core.frame.DataFrame:

    
    cooler_object : cooler.api.Cooler = cooler.Cooler(cool_file_path)
    resolution = cooler_object.binsize

    promoter_enhancer_dataframe = files.extract_pls_els_from_bed(bed_file_path)

    print("Size of promoter_enhancer_dataframe:",
            round(sys.getsizeof(promoter_enhancer_dataframe)/1048576,2),"MB", "or",round(sys.getsizeof(promoter_enhancer_dataframe),2)/1073741824,"GB")

    mcool_interaction_dataframe = note_interactions(cooler_object,promoter_enhancer_dataframe, chrom_name, start = start, end = end, resolution = resolution)

    return mcool_interaction_dataframe


def _note(cooler_file_path : str, bed_file_path :str, chrom_name :str, start : int, end : int, prev_end : int) -> pd.DataFrame:
    """ Function for use in multiprocessing in note_promoter_enhancer_interactions_threaded. 
        Finds all promoter enhancer interactions on a chromosome within a region. 

    Args:
        cooler_file_path (str): path to cooler file
        bed_file_path (str): path to bed file
        chrom_name (str): name of chromosome
        start (int): start of region to find promoter / enhancer region
        end (int): start of region to find promoter / enhancer region
        prev_end (int): if both the promoter and enhancer are positioned below this value, don't register interaction

    Returns:
        pd.Dataframe: Dataframe with promoter enhancer interactions.
    """

    print("Starting _note() process.")

    cooler_object = cooler.Cooler(cooler_file_path)
    cooler_2Dselector = cooler_object.matrix(balance=False, as_pixels=True, join=True)
    chrom_size = cooler_object.chromsizes[chrom_name]
    resolution = cooler_object.binsize

    promoter_enhancer_dataframe = files.extract_pls_els_from_bed(bed_file_path)

    new_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

    print("Fetching dataframe related to", chrom_name, "starting at", start, "ending at", end)
    cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
    print(f"Fetched dataframe for {chrom_name}:{start}-{end}. Total size: {sys.getsizeof(cooler_dataframe)/1_048_576} MB")

    ## *  Copy part of dataframe with relevant chrom_name 

    # TODO: Reduce dataframe to region within start and end variables
    reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["chrom"] == chrom_name]
    reduced_dataframe = reduced_dataframe.loc[reduced_dataframe["chromEnd"] >= start].loc[reduced_dataframe["chromStart"] <= end]
    
    ## * split dataframe into pls and els dataframes
    promoter_dataframe, enhancer_dataframe = split_df_to_pls_els(reduced_dataframe)

    ## * Iterate through promoters and find bins of correlating enhancers
    print("Finding correlating promoter and enhancer bins.")
    for row in promoter_dataframe.itertuples():
        prom_start = row[2]
        prom_end = row[3]
        prom_name = row[4]
        score = row[5] #!Unused
        strand = row[6] #!Unused
        ## * Further indexes should not be relevant

        prom_rounded_down = round_up_and_down(prom_start,resolution)[0]
        prom_rounded_up = round_up_and_down(prom_end,resolution)[1]


        ## * Getting all bins where the PLS is within the range of the first region of bin
        promoter_interaction_dataframe1 = cooler_dataframe.loc[cooler_dataframe["start1"] >= prom_rounded_down].loc[cooler_dataframe["end1"] <= prom_rounded_up]
        ## * Getting all bins where the PLS is within the range of the second region of bin
        promoter_interaction_dataframe2 = cooler_dataframe.loc[cooler_dataframe["start2"] >= prom_rounded_down].loc[cooler_dataframe["end2"] <= prom_rounded_up]


    ## * Iterating through all rows to find relevant enhancers that are in second region, while promoter is in first region
        for row2 in promoter_interaction_dataframe1.itertuples():
            enh_rounded_down = row2[5]
            enh_rounded_up = row2[6]
            count = row2[7]
                        
            ## * If we are handling promoters and enhancers handled by previous thread, continue
            if prev_end:
                if prev_end > prom_rounded_up and prev_end > enh_rounded_up :
                    continue
            
            ## * Finding all enhancers that are within bins that interacted with the PLS' bin
            enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]


            ## * Iterate through all relevant enhancers and insert into new dataframe
            for row3 in enhancer_hits_dataframe.itertuples():
                enh_start = row3[2]
                enh_end = row3[3]
                enh_name = row3[4]
                
                # * Skip if we already found interaction for this enhancer and promoter pair
                already_exists_dataframe = new_dataframe.loc[(new_dataframe['enh_name'] == enh_name) & (new_dataframe['prom_name'] == prom_name)]#f'enh_name == {enh_name} and prom_name == {prom_name}')
                if len(already_exists_dataframe):
                    continue


                ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
                diff = 0
                if prom_end < enh_start: 
                    diff = abs(enh_start - prom_end)
                    
                elif enh_end < prom_start:
                    diff = abs(prom_start - enh_end)

                if diff > DISTANCE_LIMIT: 
                    continue

                input_enh_rounded_down = round_up_and_down(enh_start,resolution)[0]
                input_enh_rounded_up = round_up_and_down(enh_end,resolution)[1]
                input_prom_rounded_down = prom_rounded_down
                input_prom_rounded_up = prom_rounded_up
                input_count = count

                # * If one of the elements span multiple bins, we take the bin with the highest count. If multiple, pick first.
                if prom_rounded_up - prom_rounded_down > resolution or input_enh_rounded_up - input_enh_rounded_down > resolution: 
                    possible_bins_dataframe = cooler_dataframe.query(f'start1 >= {prom_rounded_down} and end1 <= {prom_rounded_up} and start2 >= {enh_rounded_down} and end2 <= {enh_rounded_up}')
                    
                    if len(possible_bins_dataframe):
                        
                        highest_count = np.max(np.array(possible_bins_dataframe['count']))
                        possible_bins_dataframe = possible_bins_dataframe.loc[possible_bins_dataframe['count'] == highest_count]
                        input_prom_rounded_down = np.array(possible_bins_dataframe['start1'])[0]
                        input_prom_rounded_up = np.array(possible_bins_dataframe['end1'])[0]
                        input_enh_rounded_down = np.array(possible_bins_dataframe['start2'])[0]
                        input_enh_rounded_up = np.array(possible_bins_dataframe['end2'])[0]
                        input_count = np.array(possible_bins_dataframe['count'])[0]


                ## * Last count is -1 as we currently dont have data to insert
                internal_columns = DATAFRAME_COLUMNS_INTERNAL
                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, input_enh_rounded_down, input_enh_rounded_up, 
                            input_prom_rounded_down, input_prom_rounded_up, input_count, -1]]
                

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])

        ## * Iterating through all rows to find relevant enhancers that are in first region, while promoter is in second region        
        for row2 in promoter_interaction_dataframe2.itertuples():
            enh_rounded_down = row2[2]
            enh_rounded_up = row2[3]
            count = row2[7]

            
            ## * If we are handling promoters and enhancers handled by previous thread, go to next promoter
            if prev_end:
                if prev_end > prom_rounded_up and prev_end > enh_rounded_up and prev_end == 6_000_000:
                    continue

            ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
            enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]

            ## * Iterate through all relevant enhancers and insert into new dataframe
            for row3 in enhancer_hits_dataframe.itertuples():
                enh_start = row3[2]
                enh_end = row3[3]
                enh_name = row3[4]
                internal_columns = DATAFRAME_COLUMNS_INTERNAL
                
                # * Skip if we already found interaction for this enhancer and promoter pair (Occurs on border between bins)
                already_exists_dataframe = new_dataframe.loc[(new_dataframe['enh_name'] == enh_name) & (new_dataframe['prom_name'] == prom_name)]#f'enh_name == {enh_name} and prom_name == {prom_name}')
                if len(already_exists_dataframe):
                    continue

                ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
                diff = 0
                if prom_end < enh_start: 
                    diff = abs(enh_start - prom_end)
                    
                elif enh_end < prom_start:
                    diff = abs(prom_start - enh_end)

                if diff > DISTANCE_LIMIT: 
                    # print(f"Diff {diff} is higher than {DISTANCE_LIMIT}. Skipping.")
                    continue

                input_enh_rounded_down = round_up_and_down(enh_start,resolution)[0]
                input_enh_rounded_up = round_up_and_down(enh_end,resolution)[1]
                input_prom_rounded_down = prom_rounded_down
                input_prom_rounded_up = prom_rounded_up
                input_count = count

                # * If one of the elements span multiple bins, we take the bin with the highest count. If multiple, pick first.
                if prom_rounded_up - prom_rounded_down > resolution or input_enh_rounded_up - input_enh_rounded_down > resolution: 
                    possible_bins_dataframe = cooler_dataframe.query(f'start1 >= {enh_rounded_down} and end1 <= {enh_rounded_up} and start2 >= {prom_rounded_down} and end2 <= {prom_rounded_up}')

                    if len(possible_bins_dataframe):

                        highest_count = np.max(np.array(possible_bins_dataframe['count']))
                        possible_bins_dataframe = possible_bins_dataframe.loc[possible_bins_dataframe['count'] == highest_count]
                        input_enh_rounded_down = np.array(possible_bins_dataframe['start1'])[0]
                        input_enh_rounded_up = np.array(possible_bins_dataframe['end1'])[0]
                        input_prom_rounded_down = np.array(possible_bins_dataframe['start2'])[0]
                        input_prom_rounded_up = np.array(possible_bins_dataframe['end2'])[0]
                        input_count = np.array(possible_bins_dataframe['count'])[0]
                
                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, input_enh_rounded_down, input_enh_rounded_up, 
                            input_prom_rounded_down, input_prom_rounded_up, input_count, -1]]

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])

    print("Dropping duplicates and sorting dataframe.")
    new_dataframe = new_dataframe.drop_duplicates().sort_values(by=['prom_start','bin2_start','enh_start','bin1_start']).reset_index(drop=True)
    print(f"Returning _note() thread with params start={start},end={end},prev_end={prev_end}")
    return new_dataframe

def note_promoter_enhancer_interactions_multiprocess(cooler_file_path : str,
                                                bed_file_path : str,
                                                chrom_name : str,
                                                start : int = False,
                                                end : int = False,
                                                max_distance : int = 3_000_000,
                                                workers = 10) -> pd.core.frame.DataFrame:
    """Creates a dataframe with all promoter and enhancer interactions and their contact frequency. Includes the regulatory element name,
    bins and count. 
    TODO: Finish documentation.
    Args:
        cooler_file_path (str): _description_
        bed_file_path (str): _description_
        chrom_name (str): _description_
        start (int, optional): _description_. Defaults to False.
        end (int, optional): _description_. Defaults to False.
        max_distance (int, optional): _description_. Defaults to 3_000_000.
        workers (int, optional): _description_. Defaults to 10.

    Returns:
        pd.core.frame.DataFrame: _description_
    """
    
    dataframe_to_return : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

    cooler_object = cooler.Cooler(cooler_file_path)
    chrom_size = cooler_object.chromsizes[chrom_name]

    thread_ranges = np.empty(shape=(0,2),dtype=int)

    if end: 
        total_size = end - start
    else: 
        total_size = chrom_size
        start = 0

    covered_size : int = 0
    while covered_size < total_size:
        
        if covered_size == 0: region_start = start + covered_size

        # * We want some overlap between the regions, else we wont see contacts between them
        # * E.g. If previous region we added was 3Mbp - 6Mbp, the next one is 1.5Mbp - 7.5Mbp
        else: region_start = covered_size - max_distance
        region_end = region_start + max_distance * 2
        
        # * Don't move past the edge of the chromosome
        if region_end > total_size:
            region_end = total_size

        # * A region that a single thread will handle
        one_range = np.array([[region_start,region_end]])
        # * Collect all regions in an array
        thread_ranges = np.append(thread_ranges,one_range,axis=0)
        covered_size = region_end


    print(f'Attempting to start threads')
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:

        futures = []
        prev_end : int = False
        for one_range in thread_ranges:
            this_start = one_range[0]
            this_end = one_range[1]
            print(f"Submitting one process with function {_note.__name__} and args {this_start}, {this_end}, {prev_end}.")

            futures.append(executor.submit(_note, cooler_file_path, bed_file_path, chrom_name, this_start, this_end,
                            prev_end))
            prev_end = this_end

        print(f'All threads submitted.')
        print(f'Waiting for all threads to finish.')
        print(f'Number of processes: {len(futures)}')
        concurrent.futures.wait(futures)
        print(f'All subprocesses done.')
        
        print(f'Concatenating all dataframes.')
        for f in futures:
            returned_df = f.result()
            dataframe_to_return = pd.concat([dataframe_to_return,returned_df],ignore_index=True)
            dataframe_to_return.drop_duplicates()


        print(f'Dropping duplicates, sorting and reindexing dataframe...')
        dataframe_to_return = dataframe_to_return.drop_duplicates().sort_values(by=['prom_start','bin2_start','enh_start','bin1_start']).reset_index(drop=True)

        print("Dataframe stats:")
        print("Len:",len(dataframe_to_return))
        print(dataframe_to_return)

        return dataframe_to_return


def note_promoter_enhancer_interactions_genomewide(cooler_file_path : str,
                                                bed_file_path : str,
                                                output_path : str,
                                                cache_dataframe : bool = True, 
                                                max_distance : int = 3_000_000,
                                                workers = 10) -> pd.DataFrame:
    """Calculate promoter / enhancer interaction statistics for all 

    Args:
        cooler_file_path (str): path to .cool file.
        bed_file_path (str): path to .bed file with PLS and ELS elements. 
        output_path (str): filepath and name to save output.
        cache_dataframe (bool, optional): Wether to cache dataframes made underway or not. Can be useful if the process crashes. Defaults to True.
        max_distance (int, optional): maximum distance at which to register contacts. Defaults to 3_000_000.
        workers (int, optional): number of subprocess workers. Defaults to 10.

    Returns:
        pd.DatFrame: dataframe with data for all 
    """

    cooler_object = cooler.Cooler(cooler_file_path)
    chrom_names = cooler_object.chromnames
    
    files.create_dir_for_path(output_path)

    output_path_list = output_path.split("/")

    output_folder = ""
    if output_path_list[-1:] == "": 
        raise NameError(f"Please specify a path to a filename. {output_path} is a directory.")
        output_folder = output_path
    else:
        output_folder = "/".join(output_path_list[:-1]) + "/"
    cache_folder = output_folder + "cache/"
    files.create_dir_for_path(cache_folder)

    if os.path.isfile(output_path):
        print(f"File already exists at {output_path}. Attempting to load.")
        return files.load_dataframe(output_path)


    chrom_wide_dataframe = pd.DataFrame(columns=Constants.DATAFRAME_COLUMNS_INTERNAL)

    for name in chrom_names:
        print(f"Noting promoter and enhancer interactions for {name}.")
        
        
        
        dataframe_name = output_path_list[-1:] 
        dataframe_name = dataframe_name[0] + "." + name + ".csv" 
        dataframe_path = cache_folder + dataframe_name
        

        loaded_dataframe = files.load_dataframe(dataframe_path)
        if type(loaded_dataframe) == pd.DataFrame:
            print(f"Found cached dataframe at {dataframe_path}.")
            returned_dataframe = loaded_dataframe
        else:
            returned_dataframe = note_promoter_enhancer_interactions_multiprocess(cooler_file_path =cooler_file_path, 
                                                            bed_file_path = bed_file_path,
                                                            chrom_name= name,
                                                            max_distance = max_distance,
                                                            workers = workers)
            
            if cache_dataframe: files.save_dataframe(returned_dataframe,dataframe_path)

        chrom_wide_dataframe = pd.concat([chrom_wide_dataframe, returned_dataframe], ignore_index=True)

        
    files.save_dataframe(chrom_wide_dataframe, output_path)
    return chrom_wide_dataframe
    
def check_if_column_sorted (dataframe : pd.DataFrame,
                            col_name : str):
        column_as_numpy = np.array(dataframe[col_name])

        is_sorted = lambda a: np.all(a[:-1] <= a[1:])

        return is_sorted(column_as_numpy)





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
        if (dataframe.columns == DATAFRAME_COLUMNS_COOLER).all() and not column_set:
            columns = DATAFRAME_COLUMNS_COOLER
            bin1_start_label = "start1"
            bin1_end_label = "end1"
            bin2_start_label = "start2"
            bin2_end_label = "end2"
            count_label = "count"
            column_set = True
    except ValueError: # If column shape mismatch
        pass
    try:
        if (dataframe.columns == DATAFRAME_COLUMNS_INTERNAL).all() and not column_set:
            columns = DATAFRAME_COLUMNS_INTERNAL
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
    type_bool_df = pe_df[DATAFRAME_COLUMNS_BED[9]].str.contains(regex)
    
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

# Helper method
def split_df_to_pls_els(promoter_enhancer_dataframe : pd.core.frame.DataFrame) -> tuple:
    """Splits dataframe into two new dataframes. One containing promoter-like elements and one containing enhancer-like.
    Requires the 'type' column. 

    Args:
        promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe to split. Requires 'type' column

    Returns:
        tuple: tuple with promoter and enhancer dataframe
    """
    print("Splitting dataframe into promoter and enhancer dataframe. ")

    pls_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["type"].str.contains("PLS")]
    els_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["type"].str.contains("ELS")]
    return (pls_dataframe, els_dataframe)
