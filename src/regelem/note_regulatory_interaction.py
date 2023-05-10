import os
import warnings
import cooler                   #type:ignore
import matplotlib as mpl        #type:ignore
import matplotlib.pyplot as plt #type:ignore
import numpy as np              #type:ignore
import pandas as pd             #type:ignore
import sys 
import concurrent.futures

#from os.path import dirname, basename, isfile, join

import Constants
import files
import utils
import dataframe_functions

#! TODO: There is a LOT of room for optimisation here. Do this if you have nothing better to do (so never). 

def _note(cooler_file_path : str, bed_file_path :str, chrom_name :str, 
            start : int, end : int, prev_end : int, balanced : bool) -> pd.DataFrame:
    """ Function for use in multiprocessing. 
        Finds all promoter enhancer interactions within a region and returns a dataframe. 

    Args:
        cooler_file_path (str): path to cooler file
        bed_file_path (str): path to bed file
        chrom_name (str): name of chromosome
        start (int): start of region to find promoter / enhancer region
        end (int): start of region to find promoter / enhancer region
        prev_end (int): if both the promoter and enhancer are positioned below this value, 
                        don't register interaction. Use this if there is region overlap with 
                        a different process.

    Returns:
        pd.Dataframe: Dataframe with promoter enhancer interactions.
    """

    print("Starting _note() process.")

    cooler_object = cooler.Cooler(cooler_file_path)
    cooler_2Dselector = cooler_object.matrix(balance=balanced, as_pixels=True, join=True)
    chrom_size = cooler_object.chromsizes[chrom_name]
    resolution = cooler_object.binsize

    promoter_enhancer_dataframe = files.extract_pls_els_from_bed(bed_file_path)

    noted_interactions_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = Constants.DATAFRAME_COLUMNS_INTERNAL)

    print("Fetching dataframe related to", chrom_name, "starting at", start, "ending at", end)
    cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
    print(f"Fetched dataframe for {chrom_name}:{start}-{end}. Total size: {sys.getsizeof(cooler_dataframe)/1_048_576} MB")

    ## *  Copy part of dataframe with relevant chrom_name 

    # TODO: Reduce dataframe to region within start and end variables
    reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["chrom"] == chrom_name]
    reduced_dataframe = reduced_dataframe.loc[reduced_dataframe["chromEnd"] >= start].loc[reduced_dataframe["chromStart"] <= end]
    

    ## * split dataframe into pls and els dataframes
    promoter_dataframe, enhancer_dataframe = dataframe_functions.split_df_to_pls_els(reduced_dataframe)

    ## * Iterate through promoters and find bins of correlating enhancers
    print("Finding correlating promoter and enhancer bins.")
    for row in promoter_dataframe.itertuples():
        prom_start = row[2]
        prom_end = row[3]
        prom_name = row[4]
        ## * Further indexes should not be relevant

        ## * Find bin start and end
        prom_start_rounded_down = utils.round_up_and_down(prom_start,resolution)[0]
        prom_end_rounded_up = utils.round_up_and_down(prom_end,resolution)[1]


        ## * Getting all bins where the PLS is within the range of the first region of bin
        promoter_interaction_dataframe1 = cooler_dataframe.loc[cooler_dataframe["start1"] >= prom_start_rounded_down].loc[cooler_dataframe["end1"] <= prom_end_rounded_up]
        ## * Getting all bins where the PLS is within the range of the second region of bin
        promoter_interaction_dataframe2 = cooler_dataframe.loc[cooler_dataframe["start2"] >= prom_start_rounded_down].loc[cooler_dataframe["end2"] <= prom_end_rounded_up]

        def _iterate_and_note(dataframe : pd.DataFrame, promoter_first : bool):
            """Helper function
                Iterate through all bins in dataframe that promoter interacted with, find all bins that contains an enhancer, add to new dataframe and return.

            Args:
                dataframe (pd.DataFrame): dataframe to iterate through
                promoter_first (bool): if the promoter is in the first or second bin. 

            """

            new_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = Constants.DATAFRAME_COLUMNS_INTERNAL)


            if promoter_first:
                promoter_start_col_index = 2
                enhancer_start_col_index = 5
            else:
                promoter_start_col_index = 5
                enhancer_start_col_index = 2
            counter_index = 7
            if balanced: balanced_index = 8

            for iter_tuple in dataframe.itertuples():
                prom_bin_start = iter_tuple[promoter_start_col_index]
                prom_bin_end = iter_tuple[promoter_start_col_index + 1]
                enh_bin_start = iter_tuple[enhancer_start_col_index]
                enh_bin_end = iter_tuple[enhancer_start_col_index + 1]
                count = iter_tuple[counter_index]
                if balanced: balanced_val = iter_tuple[balanced_index]
                else: balanced_val = -1

                ## * If we are handling promoters and enhancers handled by previous process, continue
                if prev_end:
                    if prev_end > prom_end_rounded_up and prev_end > enh_bin_end :
                        continue
                
                ## * Finding all enhancers that are within bins that interacted with the PLS bin
                enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_bin_end].loc[enhancer_dataframe["chromEnd"] > enh_bin_start]

                for enhancer_tuple in enhancer_hits_dataframe.itertuples():
                    enh_start = enhancer_tuple[2]
                    enh_end = enhancer_tuple[3]
                    enh_name = enhancer_tuple[4]


                    ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
                    diff = 0
                    if prom_end < enh_start: 
                        diff = abs(enh_start - prom_end)
                        
                    elif enh_end < prom_start:
                        diff = abs(prom_start - enh_end)

                    if diff > Constants.DISTANCE_LIMIT: 
                        continue

                    input_list = [[chrom_name, 
                                    enh_start, enh_end, 
                                    prom_start, prom_end,
                                    enh_name, prom_name, 
                                    enh_bin_start, enh_bin_end,
                                    prom_bin_start, prom_bin_end,
                                    count, balanced_val]]

                    input_df = pd.DataFrame(input_list, columns = Constants.DATAFRAME_COLUMNS_INTERNAL)
                    new_dataframe = pd.concat([new_dataframe,input_df],ignore_index=True)
            return new_dataframe

        noted_interactions_dataframe = pd.concat([noted_interactions_dataframe,_iterate_and_note(promoter_interaction_dataframe1, True)])
        noted_interactions_dataframe = pd.concat([noted_interactions_dataframe,_iterate_and_note(promoter_interaction_dataframe2, False)])
        

    print("Dropping duplicates and sorting dataframe.")
    noted_interactions_dataframe = noted_interactions_dataframe.drop_duplicates().sort_values(by=['prom_start','bin2_start','enh_start','bin1_start']).reset_index(drop=True)
    print(f"Returning _note() thread with params start={start},end={end},prev_end={prev_end}")


    return noted_interactions_dataframe


def note_prom_enh_inter_multiprocess(cooler_file_path : str,
                                                bed_file_path : str,
                                                chrom_name : str,
                                                start : int = False,
                                                end : int = False,
                                                max_distance : int = 3_000_000,
                                                workers = 10,
                                                balanced : bool = True) -> pd.core.frame.DataFrame:
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
    
    dataframe_to_return : pd.core.frame.DataFrame = pd.DataFrame(columns = Constants.DATAFRAME_COLUMNS_INTERNAL)

    cooler_object = cooler.Cooler(cooler_file_path)
    chrom_size = cooler_object.chromsizes[chrom_name]

    region_ranges = np.empty(shape=(0,2),dtype=int)

    if start and end: 
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
        region_ranges = np.append(region_ranges,one_range,axis=0)
        covered_size = region_end


    print(f'Attempting to start subprocesses')
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:

        futures = []
        prev_end : int = False
        for one_range in region_ranges:
            this_start = one_range[0]
            this_end = one_range[1]
            print(f"Submitting one process with function {_note.__name__} and args {this_start}, {this_end}, {prev_end}.")

            futures.append(executor.submit(_note, cooler_file_path, bed_file_path, chrom_name, this_start, this_end,
                            prev_end, balanced))
            prev_end = this_end

        print(f'All threads submitted.')
        print(f'Waiting for all subprocesses to finish.')
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


def note_prom_enh_inter_gwide(cooler_file_path : str,
                                                bed_file_path : str,
                                                output_path : str,
                                                balanced : bool = True,
                                                cache_dataframe : bool = True, 
                                                max_distance : int = 3_000_000,
                                                workers = 10) -> pd.DataFrame:
    """Find all promoter-enhancer interactions for all chromosomes listed in a cooler file based on regulatory elements of a bed file.

    Args:
        cooler_file_path (str): path to .cool file.
        bed_file_path (str): path to .bed file with PLS and ELS elements. 
        output_path (str): filepath to save output.
        cache_dataframe (bool, optional): If true, cache smaller dataframes while processing. Can save time if the process is interrupted. Defaults to True.
        max_distance (int, optional): maximum distance to register contacts. Defaults to 3_000_000.
        workers (int, optional): number of subprocess workers to be used. Defaults to 10.
        balanced (bool, optional): If the registered counts should be balanced or not. Defaults to True.

    Raises:
        ValueError: Incorrect output_path value. 

    Returns:
        pd.DataFrame: Dataframe containing all promoter-enhancer interactions.
    """
    
    if os.path.isfile(output_path):
        print(f"File already exists at {output_path}. Attempting to load.")
        return files.load_dataframe(output_path)

    cooler_object = cooler.Cooler(cooler_file_path)
    chrom_names = cooler_object.chromnames
    
    files.create_dir_for_path(output_path)

    output_path_list = output_path.split("/")

    output_folder = ""
    if output_path_list[-1:] == "": 
        raise ValueError(f"Please specify a path to a filename. {output_path} is a directory.")
        output_folder = output_path
    else:
        output_folder = "/".join(output_path_list[:-1]) + "/"
    cache_folder = output_folder + "cache/"
    files.create_dir_for_path(cache_folder)



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
            returned_dataframe = note_prom_enh_inter_multiprocess(cooler_file_path =cooler_file_path, 
                                                            bed_file_path = bed_file_path,
                                                            chrom_name= name,
                                                            max_distance = max_distance,
                                                            workers = workers,
                                                            balanced = balanced)
            
            if cache_dataframe: files.save_dataframe(returned_dataframe,dataframe_path)

        chrom_wide_dataframe = pd.concat([chrom_wide_dataframe, returned_dataframe], ignore_index=True)

    chrom_wide_dataframe = chrom_wide_dataframe.drop_duplicates().sort_values(by=['prom_start','bin2_start','enh_start','bin1_start']).reset_index(drop=True)
    files.save_dataframe(chrom_wide_dataframe, output_path)
    return chrom_wide_dataframe








