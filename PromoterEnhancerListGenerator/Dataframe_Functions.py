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
from . import (File_Functions, Functions, General_Functions, Constants)
from Constants import *
from .File_Functions import *
from .Functions import *
from .General_Functions import *


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
def create_regElem_df_of_cool(
                        cool_file_path : str,
                        bed_file_path : str,
                        chrom_name : str, 
                        resolution : int = False,
                        start : int = None,
                        end : int = None,
                        max_memory : int = None) -> pd.core.frame.DataFrame:
    
    cooler_object : cooler.api.Cooler = cooler.Cooler(cool_file_path)
    resolution = cooler_object.binsize

    promoter_enhancer_dataframe = extract_pls_els_from_bed(bed_file_path)

    print("Size of promoter_enhancer_dataframe:",
            round(sys.getsizeof(promoter_enhancer_dataframe)/1048576,2),"MB", "or",round(sys.getsizeof(promoter_enhancer_dataframe),2)/1073741824,"GB")

    mcool_interaction_dataframe = note_interactions(cooler_object,promoter_enhancer_dataframe, chrom_name, start = start, end = end, resolution = resolution)

    return mcool_interaction_dataframe

def note_promoter_enhancer_interactions_threaded(cooler_file_path : str,
                                                promoter_enhancer_dataframe : pd.core.frame.DataFrame,
                                                chrom_name : str,
                                                start : int = False,
                                                end : int = False,
                                                max_distance : int = 3_000_000,
                                                workers = 10) -> pd.core.frame.DataFrame:
    
    collection_lock = threading.Lock()

    def _collect_threads_to_dataframe(dataframe_to_add : pd.DataFrame):
        with collection_lock:
            dataframe_to_return = dataframe_to_return.concat([dataframe_to_return,dataframe_to_add],ignore_index=True)
            dataframe_to_return.drop_duplicates()
    

    def _note(start : int, end : int, prev_end : int = False):

        print(f"Starting _note() thread with params start={start},end={end},prev_end={prev_end}")

        new_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

        print("Fetching dataframe related to", chrom_name, "starting at", start, "ending at", end)
        cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
        print(f"Fetched dataframe for {chrom_name}:{start}-{end}. Total size: {sys.getsizeof(cooler_dataframe)/1_048_576}MB")

        ## *  Copy part of dataframe with relevant chrom_name 
        reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[DATAFRAME_COLUMNS_BED[0]] == chrom_name]
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

            prom_start_rounded_down = round_up_and_down(prom_start,resolution)[0]
            prom_end_rounded_up = round_up_and_down(prom_end,resolution)[1]
            ## * Getting all rows of interactions where the PLS is within the bin range
            promoter_interaction_dataframe = cooler_dataframe.loc[cooler_dataframe["start1"] == prom_start_rounded_down]
            if prom_end_rounded_up - resolution != prom_start_rounded_down: 
                rest_of_promoter_interaction_dataframe = cooler_dataframe.loc[cooler_dataframe["end1"] == prom_end_rounded_up]
                promoter_interaction_dataframe = pd.concat([promoter_interaction_dataframe,rest_of_promoter_interaction_dataframe],
                                                            ignore_index=True,
                                                            )
                duplicates = promoter_interaction_dataframe.duplicated(keep=False)
                if duplicates.any():
                    print(promoter_interaction_dataframe[duplicates])
                    exit()

            ## * Iterating through all rows to find relevant enhancers
            for row2 in promoter_interaction_dataframe.itertuples():
                enh_rounded_down = row2[5]
                enh_rounded_up = row2[6]
                modle_count = row2[7]
                
                ## * If distance between promoter/enhancer is longer than max_distance, go to next promoter
                diff = abs(prom_start_rounded_down - enh_rounded_down)
                if diff > max_distance: break
                
                ## * If we are handling promoters and enhancers handled by previous thread, go to next promoter
                if prev_end:
                    if prev_end > prom_end and prev_end > enh_end:
                        break
                
                ## * Finding all enhancers that are within bins that interacted with the PLS' bin
                enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] >= enh_rounded_down].loc[enhancer_dataframe["chromStart"] <= enh_rounded_up]

                ## * Iterate through all relevant enhancers and insert into new dataframe
                for row3 in enhancer_hits_dataframe.itertuples():
                    enh_start = row3[2]
                    enh_end = row3[3]
                    enh_name = row3[4]
                    
                    ## * Last count is -1 as we currently dont have data to insert
                    internal_columns = DATAFRAME_COLUMNS_INTERNAL

                    ## * Make a list out of data to insert into new dataframe
                    ## * has to be a list of list
                    input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                                enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                                prom_start_rounded_down, prom_end_rounded_up, modle_count, -1]]

                    input_df = pd.DataFrame(input_list, columns = internal_columns)
                    new_dataframe = pd.concat([new_dataframe,input_df])


        print("Sorting dataframe and dropping duplicates.")
        new_dataframe = new_dataframe.drop_duplicates()
        _collect_threads_to_dataframe(new_dataframe)
        
    
    dataframe_to_return : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

    cooler_object = cooler.Cooler(cooler_file_path)
    cooler_2Dselector = cooler_object.matrix(balance=False, as_pixels=True, join=True)
    chrom_size = cooler_object.chromsizes[chrom_name]
    resolution = cooler_object.binsize

    thread_ranges = np.empty(shape=(0,2),dtype=int)

    if start and end: 
        total_size = end - start
    else: 
        total_size = chrom_size

    covered_size : int = 0
    while covered_size < total_size:
        
        
        if covered_size == 0: region_start = start + covered_size
        # * We want some overlap between the regions, else we wont see contacts between them
        else: region_start = covered_size - max_distance

        region_end = region_start + max_distance * 2
        
        if region_end > total_size:
            region_end = total_size


        one_range = np.array([[region_start,region_end]])
        thread_ranges = np.append(thread_ranges,one_range,axis=0)
        covered_size = region_end


    print(f'Attempting ')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:

        futures = []
        prev_end : int = False
        for one_range in thread_ranges:
            this_start = one_range[0]
            this_end = one_range[1]

            print(f"Submitting one thread with function {_note.__name__} and args {this_start}, {this_end}, {prev_end}.")
            futures.append(executor.submit(_note,this_start,this_end,prev_end))
            prev_end = this_end

        print(f'All threads submitted.')
        print(f'Waiting for all threads to finish.')
        print(f'len{futures}')
        concurrent.futures.wait(futures)
        
        for f in futures:
            print(f.result())

        exit()
        print(f'All threads done.')
        print(f'Dropping duplicates, sorting and reindexing dataframe...')
        dataframe_to_return.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

        print("Dataframe stats:")
        print("Len:",len(dataframe_to_return))
        print(dataframe_to_return)

        return dataframe_to_return

    


def note_interactions(  cool_file, 
                        promoter_enhancer_dataframe : pd.core.frame.DataFrame, 
                        chrom_name : str,
                        start : int = None,
                        end : int = None, 
                        resolution : int = False,
                        max_distance : int = False) -> pd.core.frame.DataFrame:
    """Looks through pixels in a cooler dataframe of a given resolution and finds those who include interactions
    between a promoter-like element and a enhancer-like element. Returns a dataframe. 


    Args:
        cooler_2Dselector (cooler.core.RangeSelector2D): a cooler matrix selector
        promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe with promoter/enhancer regions
        chrom_name (str): chromosome we want to look at
        start (int, optional): start of frame we want to look at. Defaults to None.
        end (int, optional): end of frame we want to look at. Defaults to None.
        resolution (int, optional): Deprecated.
        #TODO Make resolution dynamic like in note_validation
    Returns:
        pd.core.frame.DataFrame: Dataframe with counts of proposed promoter/enhancer interactions. Columns = DATAFRAME_COLUMNS_INTERNAL
    """

    print("\n\nNoting interactions for chrom:", chrom_name, "at resolution", resolution, "...")

    new_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

    cooler_2Dselector = cool_file.matrix(balance=False, as_pixels=True, join=True)
    chrom_size = cool_file.chromsizes[chrom_name]
    resolution = cool_file.binsize

    # TODO Make it so we can fetch full chromosome without running out of memory.


    # TODO fetch numpy array and index manually. Roberto says this is faster than fetching dataframe
    if start and end:
        print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
        cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
    else:
        print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
        print("Note: This may eat a lot of memory.")
        cooler_dataframe = cooler_2Dselector.fetch((chrom_name))

    ## *  Copy part of dataframe with relevant chrom_name 
    reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[DATAFRAME_COLUMNS_BED[0]] == chrom_name]
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
    
        prom_rounded_down, prom_rounded_up = round_up_and_down(prom_start,resolution)
        ## * Getting all interaction rows where the PLS is within the bin range
        promoter_interaction_dataframe = cooler_dataframe.loc[cooler_dataframe["start1"] == prom_rounded_down]
        ## * Iterating through all rows to find relevant enhancers
        for row2 in promoter_interaction_dataframe.itertuples():
            enh_rounded_down = row2[5]
            enh_rounded_up = row2[6]
            modle_count = row2[7]
            
            ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
            diff = abs(prom_rounded_down - enh_rounded_down)
            if diff > DISTANCE_LIMIT: continue

            ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
            enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] >= enh_rounded_down].loc[enhancer_dataframe["chromStart"] <= enh_rounded_up]

            ## * Iterate through all relevant enhancers and insert into new dataframe
            for row3 in enhancer_hits_dataframe.itertuples():
                enh_start = row3[2]
                enh_end = row3[3]
                enh_name = row3[4]
                
                ## * Last count is -1 as we currently dont have data to insert
                internal_columns = DATAFRAME_COLUMNS_INTERNAL

                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                            prom_rounded_down, prom_rounded_up, modle_count, -1]]

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])


    print("Sorting dataframe and dropping duplicates.")
    new_dataframe = new_dataframe.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

    return new_dataframe

def note_validation_data(cooler_2Dselector : cooler.core.RangeSelector2D, 
                        promoter_enhancer_dataframe: pd.core.frame.DataFrame, 
                        modle_dataframe: pd.core.frame.DataFrame, 
                        chrom_name : str, 
                        start : int = None, 
                        end : int = None,
                        count_column : int = 1) -> pd.core.frame.DataFrame:
    """Notes promoter/enhancer interaction for a validation dataframe and adds its values to the modle dataframe.

    Args:
        cooler_2Dselector (cooler.core.RangeSelector2D): a cooler matrix selector
        promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe with promoter/enhancer regions
        modle_dataframe (pd.core.frame.DataFrame): dataframe based on data generated by modle
        chrom_name (str): chromosome we want to look at
        start (int, optional): start of frame we want to look at. Defaults to None.
        end (int, optional): end of frame we want to look at. Defaults to None.

    Returns:
        pd.core.frame.DataFrame: Modle dataframe with respective validation numbers. 
    """
    
    if start and end:
        print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
        valid_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
    else:
        print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
        print("        Note: This may eat a lot memory")
        valid_dataframe = cooler_2Dselector.fetch((chrom_name))

    resolution = valid_dataframe["start1"][0] - valid_dataframe["end1"][0]

    valid_dataframe = note_interactions(cooler_2Dselector, promoter_enhancer_dataframe, chrom_name, start, end, resolution)
    
    # Swap columns
    modle_count_column = valid_dataframe['modle_count']
    valid_count_column = valid_dataframe['valid_count']
    valid_dataframe['modle_count'] = valid_count_column
    valid_dataframe['valid_count'] = modle_count_column
    

    print("Counting up interactions from validation data.")
    # Collect all counts within 1000 res dataframe to a 5000 one
    for row in modle_dataframe.itertuples():
        row_index = row[0]
        bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
        bin1_start = row[bin1_start_index + 1]
        bin1_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_end")
        bin1_end = row[bin1_end_index + 1]
        
        bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
        bin2_start = row[bin2_start_index+ 1]
        bin2_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_end")
        bin2_end = row[bin2_end_index + 1]


        start = round_up_and_down(bin1_start,5000)


        round_towards_bin1_start = lambda param : round_up_and_down(param,5000)[0] == bin1_start
        round_towards_bin2_start = lambda param : round_up_and_down(param,5000)[0] == bin2_start

        bin2_in_range = valid_dataframe['bin2_start'].map(round_towards_bin2_start)
        filtered_dataframe = valid_dataframe.loc[bin2_in_range]

        #print(filtered_dataframe)
        bin1_in_range = filtered_dataframe['bin1_start'].map(round_towards_bin1_start)
        filtered_dataframe = filtered_dataframe.loc[bin1_in_range]

        sum_of_pixel = filtered_dataframe['valid_count'].sum()
        
        modle_dataframe.at[row_index, 'valid_count'] = sum_of_pixel


    print("Validation data added to dataframe")
    return modle_dataframe

def collect_dataframe_to_match_resolution(dataframe: pd.core.frame.DataFrame,
                                        target_resolution : int = 5000):

    ## ! Tacky, bad solution
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




def create_list_of_chrom(
                        cooler_file_path : str, 
                        bed_file_path : str, 
                        chrom_name : str, 
                        start : int = None, 
                        end : int = None, 
                        cooler_valid_file_path : str = None) -> pd.core.frame.DataFrame:




    cooler_object, cooler_2Dselector = read_cool_file(cooler_file_path)

    chrom_names : list = cooler_object.chromnames

    if chrom_name not in chrom_names:
        #TODO Raise error
        pass
    
    bed_df : pd.core.frame.DataFrame = extract_pls_els_from_bed(bed_file_path)
    promoter_enhancer_dataframe : pd.core.frame.DataFrame = filter_type_in_dataframe(bed_df)
    
    interactions_dataframe : pd.core.frame.DataFrame = note_interactions(cooler_2Dselector, 
                                                                        promoter_enhancer_dataframe, 
                                                                        chrom_name, 
                                                                        start = start, 
                                                                        end = end)
    

    if cooler_valid_file_path:
        cooler_valid_file_path_5000 = cooler_valid_file_path + "::/resolutions/5000"
        cooler_valid_file_path_1000 = cooler_valid_file_path + "::/resolutions/1000"

        cool_valid_object_5000 : cooler.api.Cooler = read_mcool_file(cooler_valid_file_path_5000)
        valid_2D_selector_5000 = cool_valid_object_5000.matrix(balance=False, as_pixels=True, join=True)

        cool_valid_object_1000 : cooler.api.Cooler = read_mcool_file(cooler_valid_file_path_1000)
        valid_2D_selector_1000 = cool_valid_object_1000.matrix(balance=False, as_pixels=True, join=True)

        interactions_dataframe = note_validation_data(valid_2D_selector_1000, promoter_enhancer_dataframe, interactions_dataframe, chrom_name, start, end)

    write_to_file_bedpe(interactions_dataframe, "defaultOutputName.txt")

    return interactions_dataframe


def create_distance_dataframe(dataframe: pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    """Creates a dataframe where each row has the distance between two bins, and the contact frequency between them"""

    print("Creating distance dataframe.")

    format1 = False
    format2 = False

    def _format1():
        distance_matrix = pd.DataFrame(columns=["distance","count"])
        for row in dataframe.itertuples():
            bin1_start_index = dataframe.columns.get_loc("start1")
            bin1_start = row[bin1_start_index + 1]
            bin2_start_index = dataframe.columns.get_loc("start2")
            bin2_start = row[bin2_start_index + 1]

            distance = abs(bin2_start - bin1_start)
            count_index = dataframe.columns.get_loc("count")
            count = row[count_index + 1]

            distance_matrix.loc[len(distance_matrix)] = [distance,count]

        distance_matrix = distance_matrix.convert_dtypes()
        return distance_matrix

    def _format2():
        distance_matrix = pd.DataFrame(columns=["distance","count"])
        for row in dataframe.itertuples():
            bin1_start_index = dataframe.columns.get_loc("prom_start")
            bin1_start = row[bin1_start_index + 1]
            bin2_start_index = dataframe.columns.get_loc("enh_start")
            bin2_start = row[bin2_start_index + 1]

            distance = abs(bin2_start - bin1_start)
            count_index = dataframe.columns.get_loc("modle_count")
            count = row[count_index + 1]

            distance_matrix.loc[len(distance_matrix)] = [distance,count]

        distance_matrix = distance_matrix.convert_dtypes()
        return distance_matrix

    if ("start1" in dataframe.columns) and ("start2" in dataframe.columns) and ("count" in dataframe.columns):
        format1 = True
    elif ("bin1_start" in dataframe.columns) and ("bin2_start" in dataframe.columns) and ("modle_count" in dataframe.columns):
        format2 = True
    else:    
        print("Invalid format.")
    
    if format1: 
        return _format1()
    if format2:
        return _format2()


def expand_distance_dataframe(dataframe : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    """ Looks for 'distance" and 'count' columsn in the given dataframe. 
    If present, expands the dataframe to also include the foolowing columns:
            'distance', 'mincount', 'maxcount,'
            'averagecount', 'mediancount', 'standarddeviation',
            'totalcount', 'numberofcounts', 'allcounts'


    """

    print("Expanding distance dataframe with more statistics.")

    # Check if dataframe already has correct format. If not, reformat. 
    if not "distance" in dataframe.columns and not "count" in dataframe.columns:
        print("NO__")

    
    # Find all unique distances
    uqdistances = np.unique(dataframe['distance'].to_numpy())
    emptyarray = [None for i in uqdistances]
    emptyarrayarray = [[] for i in uqdistances]

    length = len(uqdistances)
    # Temporary dictionary to easily make a dataframe (because pandas is a pain)
    dictionary = {'distance': uqdistances, 'mincount': emptyarray, 'maxcount': emptyarray,
            'averagecount': emptyarray, 'mediancount': emptyarray, 'standarddeviation': emptyarray,
            'totalcount': emptyarray, 'numberofcounts': emptyarray, 'allcounts': emptyarrayarray}


    # Create a dataframe with distance as the index
    new_dataframe = pd.DataFrame.from_dict(dictionary)
    new_dataframe['distance'] = new_dataframe['distance'].astype('int64')
    new_dataframe = new_dataframe.set_index('distance')

    # Iterate through the old dataframe
    for row in dataframe.itertuples():
        distance = int(row[1])
        count = row[2]

        # Insert all counts with the same distance into an array at the same index 
        new_dataframe.at[distance,'allcounts'].append(count)
    # Iterate through new dataframe and update all data based on the counts the 'allcounts' array
    for row in new_dataframe.itertuples():
        allcounts = row[DATAFRAME_COLUMNS_STATISTICAL.index('allcounts')]
        mincount = min(allcounts)
        maxcounts = max(allcounts)
        averagecount = np.mean(allcounts)
        mediancount = np.median(allcounts)
        standarddeviation = np.std(allcounts)
        totalcount = np.sum(allcounts)
        numberofcounts = len(allcounts)

        # Disables warning from the following line
        np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 
        new_dataframe.loc[row[0]] = [mincount,maxcounts,averagecount,mediancount,standarddeviation,totalcount,numberofcounts,allcounts]

    return new_dataframe

# def upscale_distance_dataframe(dataframe : pd.core.frame.DataFrame,
#                                 target_res) -> pd.core.frame.Dataframe:

                                

#         resolution = dataframe["bin1_start"][0] - dataframe["bin1_end"][0]

#     if target_resolution % resolution != 0:
#         print("target_resolution not divisible with original resolution")
#         return False

#     print("Counting up interactions from validation data")

#     dataframe_with_target_res = pd.core.frame.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)


#     ## * Collect all counts within 1000 res dataframe to a 5000 one
#     for row in dataframe.itertuples():
#         row_index = row[0]
#         bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
#         bin1_start = row[bin1_start_index + 1]
#         bin1_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_end")
#         bin1_end = row[bin1_end_index + 1]
        
#         bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
#         bin2_start = row[bin2_start_index+ 1]
#         bin2_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_end")
#         bin2_end = row[bin2_end_index + 1]

#         bin1_start_rounded, bin1_end_rounded = round_up_and_down(bin1_start,5000)
#         bin2_start_rounded, bin2_end_rounded = round_up_and_down(bin2_start,5000)


#         count_index = DATAFRAME_COLUMNS_INTERNAL.index('modle_count')
#         count = row[count_index + 1]

#         target_row = dataframe_with_target_res.loc[(dataframe_with_target_res['bin1_start']==bin1_start_rounded) & (dataframe_with_target_res['bin2_start']==bin2_start_rounded)]

#         if target_row.empty:
#                         row_list = list(row)
#                         row_list[bin1_start_index + 1] = bin1_start_rounded
#                         row_list[bin1_end_index + 1] = bin1_end_rounded
#                         row_list[bin2_start_index + 1] = bin2_start_rounded
#                         row_list[bin2_end_index + 1] = bin2_end_rounded
                        
#                         dataframe_with_target_res = pd.concat([dataframe_with_target_res, pd.core.frame.DataFrame([row_list[1:]], columns=DATAFRAME_COLUMNS_INTERNAL)]).reset_index(drop=True)
#         else:
#             indexes = dataframe_with_target_res.index[(dataframe_with_target_res['bin1_start']==bin1_start_rounded) & (dataframe_with_target_res['bin2_start']==bin2_start_rounded)]
#             if len(indexes) > 1: 
#                 print("Duplicate rows")
#                 exit()
#             index = indexes[0]
#             #if count > 100: 
#                 #print(count)
#             #print(dataframe_with_target_res.iloc[index]['modle_count'])
#             dataframe_with_target_res.iloc[index]['modle_count'] += count

#     dataframe_with_target_res = dataframe_with_target_res.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

#     return dataframe_with_target_res


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

    Args:
        promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe to split

    Returns:
        tuple: tuple with promoter and enhancer dataframe
    """
    print("Splitting dataframe into promoter and enhancer dataframe. ")

    type_index = DATAFRAME_COLUMNS_BED.index("type")
    pls_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[DATAFRAME_COLUMNS_BED[type_index]].str.contains("PLS")]
    els_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe[DATAFRAME_COLUMNS_BED[type_index]].str.contains("ELS")]
    return (pls_dataframe, els_dataframe)
