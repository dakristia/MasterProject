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
from Constants import *
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

    print(cooler_file_path)

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
            
            ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
            ## * In theory this shouldn't really happen, as most if not all promoter/enhancer interactions are within 3 Mbp of eachother. 
            diff = 0
            if prom_rounded_up < enh_rounded_down: 
                diff = abs(prom_rounded_up - enh_rounded_down)
                
            elif enh_rounded_up < prom_rounded_down:
                diff = abs(enh_rounded_up - prom_rounded_down)

            if diff > DISTANCE_LIMIT: 
                print(f"Diff {diff} is higher than {DISTANCE_LIMIT}. Skipping row.")
                continue

            
            ## * If we are handling promoters and enhancers handled by previous thread, go to next promoter
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
                
                ## * Last count is -1 as we currently dont have data to insert
                internal_columns = DATAFRAME_COLUMNS_INTERNAL

                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                            prom_rounded_down, prom_rounded_up, count, -1]]

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
                
                ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
                ## * In theory this shouldn't really happen, as most if not all promoter/enhancer interactions are within 3 Mbp of eachother. 
                diff = 0
                if prom_end < enh_start: 
                    diff = abs(enh_start - prom_end)
                    
                elif enh_end < prom_start:
                    diff = abs(prom_start - enh_end)

                if diff > DISTANCE_LIMIT: 
                    print(f"Diff {diff} is higher than {DISTANCE_LIMIT}. Skipping.")
                    continue

                ## * If any of the elements span multiple bins, we only count the bin containing the start of the bin
                if enh_rounded_up - enh_rounded_down > resolution:
                    enh_rounded_up = enh_rounded_down + resolution
                if prom_rounded_up - prom_rounded_down > resolution:
                    prom_rounded_up = prom_rounded_down + resolution

                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                            prom_rounded_down, prom_rounded_up, count, -1]]

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])

    print("Dropping duplicates and sorting dataframe.")
    new_dataframe = new_dataframe.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)
    print(f"Returning _note() thread with params start={start},end={end},prev_end={prev_end}")
    return new_dataframe



def note_promoter_enhancer_interactions_multiprocess(cooler_file_path : str,
                                                bed_file_path : str,
                                                chrom_name : str,
                                                start : int = False,
                                                end : int = False,
                                                max_distance : int = 3_000_000,
                                                workers = 10) -> pd.core.frame.DataFrame:
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
            print(f"Submitting one thread with function {_note.__name__} and args {this_start}, {this_end}, {prev_end}.")
            #futures.append(executor.submit(_note,cooler_2Dselector, promoter_enhancer_dataframe, 
            #                chrom_name, resolution, this_start,this_end,prev_end))
            futures.append(executor.submit(_note, cooler_file_path, bed_file_path, chrom_name, this_start, this_end,
                            prev_end))
            prev_end = this_end

        print(f'All threads submitted.')
        print(f'Waiting for all threads to finish.')
        print(f'Number of threads: {len(futures)}')
        concurrent.futures.wait(futures)
        print(f'All threads done.')
        
        print(f'Concatenating all dataframes.')
        for f in futures:
            returned_df = f.result()
            dataframe_to_return = pd.concat([dataframe_to_return,returned_df],ignore_index=True)
            dataframe_to_return.drop_duplicates()


        print(f'Dropping duplicates, sorting and reindexing dataframe...')
        dataframe_to_return = dataframe_to_return.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

        print("Dataframe stats:")
        print("Len:",len(dataframe_to_return))
        print(dataframe_to_return)

        return dataframe_to_return

    
def check_if_column_sorted (dataframe : pd.DataFrame,
                            col_name : str):
        column_as_numpy = np.array(dataframe[col_name])

        is_sorted = lambda a: np.all(a[:-1] <= a[1:])

        return is_sorted(column_as_numpy)



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

    # TODO fetch numpy array and index manually. Roberto says this is faster than fetching dataframe
    if end:
        print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
        cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
    else:
        print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
        print("Note: This may eat a lot of memory.")
        cooler_dataframe = cooler_2Dselector.fetch((chrom_name))

    ## *  Copy part of dataframe with relevant chrom_name 
    reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["chrom"] == chrom_name]
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

        # Get the lower end of bins    
        prom_rounded_down = round_up_and_down(prom_start,resolution)[0]
        # Get the higher end of bins
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
            

            ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
            ## * In theory this shouldn't really happen, as most if not all promoter/enhancer interactions are within 3 Mbp of eachother. 
            diff = 0
            if prom_rounded_up < enh_rounded_down: 
                diff = abs(prom_rounded_up - enh_rounded_down)
                
            elif enh_rounded_up < prom_rounded_down:
                diff = abs(enh_rounded_up - prom_rounded_down)

            if diff > DISTANCE_LIMIT: 
                print(f"Diff {diff} is higher than {DISTANCE_LIMIT}. Skipping row.")
                continue

            ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
            #enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] >= enh_rounded_down].loc[enhancer_dataframe["chromStart"] < enh_rounded_up]
            enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]

            ## * Iterate through all relevant enhancers and insert into new dataframe
            for row3 in enhancer_hits_dataframe.itertuples():
                enh_start = row3[2]
                enh_end = row3[3]
                enh_name = row3[4]
                
                ## * Last count is -1 as we currently dont have data to insert
                internal_columns = DATAFRAME_COLUMNS_INTERNAL

                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                            prom_rounded_down, prom_rounded_up, count, -1]]

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])
        
        ## * Iterating through all rows to find relevant enhancers that are in first region, while promoter is in second region        
        for row2 in promoter_interaction_dataframe2.itertuples():
            enh_rounded_down = row2[2]
            enh_rounded_up = row2[3]
            count = row2[7]



            ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
            diff = abs(prom_rounded_down - enh_rounded_down)
            if diff > DISTANCE_LIMIT: continue

            ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
            enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]

            ## * Iterate through all relevant enhancers and insert into new dataframe
            for row3 in enhancer_hits_dataframe.itertuples():
                enh_start = row3[2]
                enh_end = row3[3]
                enh_name = row3[4]

                internal_columns = DATAFRAME_COLUMNS_INTERNAL

                ## * Make a list out of data to insert into new dataframe
                ## * has to be a list of list
                #TODO: Change column order so promoters are first, and enhancers second.
                input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
                            enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
                            prom_rounded_down, prom_rounded_up, count, -1]]

                input_df = pd.DataFrame(input_list, columns = internal_columns)
                new_dataframe = pd.concat([new_dataframe,input_df])


    print("Sorting dataframe and dropping duplicates.")
    new_dataframe = new_dataframe.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

    print(new_dataframe)
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
