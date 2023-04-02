import itertools
import sys
import os
import traceback
from os.path import dirname, basename, isfile, join

#sys.path.insert(0, '..')
from Constants import *
import numpy as np
import pandas as pd
import PromoterEnhancerListGenerator as pelg
import Plotter as plotter 
import matplotlib.pyplot as plt #type:ignore
import cooler
import time
import math
from generate_plots import *
from datetime import datetime
from files import *
import threading
from threading import Thread
import concurrent.futures


input_folder = "../../input/"
output_folder = "../../output/"
predicted_output_folder = f"{output_folder}predicted/"

temp_files = "./temp/"    

def predict_average_promoter_enhancer_count_for_all_chroms_and_resolutions():
    """Naively predicts the average promoter and enhancer count for all chromosomes and resolutions. Does this by reading a bed file with PLS and ELS elements, 
    and counting any bin with a possible PLS and ELS interactions as a definit interaction.
    """
    

    max_workers = 20
    chrom_sizes_file_path = input_folder + "hg38.chrom.sizes"
    bed_file_path = input_folder + "H1-hESC.7group.bed"

    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])

    # * Get promoters and enhancers from bed file
    promoter_enhancer_dataframe = pd.read_csv(bed_file_path,delim_whitespace=True,header=None,names=DATAFRAME_COLUMNS_BED)
    # * Filter 'type' column to only include PLS, pELS or dELS (makes it easier to check this column later)
    promoter_enhancer_dataframe = pelg.filter_type_in_dataframe(promoter_enhancer_dataframe)
    # * Split dataframe into one for promoters and one for enhancers
    promoter_dataframe, enhancer_dataframe = pelg.split_df_to_pls_els(promoter_enhancer_dataframe)

    # * Array with all resolutions we want to calculate for
    resolutions = np.array([50,100,500,1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000])


    print(f"Starting {max_workers} processes with function {calculate_promoter_enhancer_bins_multiprocess.__name__}")
    #with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = []
    for row in chrom_sizes_dataframe.iloc[::-1].itertuples():
        chrom_name = row[1]
        chrom_size = row[2]
        total_parts = len(promoter_dataframe.loc[promoter_dataframe['chrom'] == chrom_name])
        for resolution in resolutions:
            out_file_name = f'{predicted_output_folder}predicted_chromosome_data_{chrom_name}_{resolution}.csv'
            if os.path.isfile(out_file_name):
                print(f'File at path {out_file_name} already exists. Skipping.')
                continue
            print(f'File {out_file_name} does NOT exist. Predict!')
            calculate_promoter_enhancer_bins_multiprocess(bed_file_path,chrom_name,chrom_size,resolution,10,total_parts,out_file_name)
        break #! Remove

    print("Done.")


def _predict(bed_file_path : str, chrom_name : str, resolution : int, total_parts : int, part : int, max_distance : int = 3_000_000):
    """Calculate possible promoter enhancer interactions based on data from a bed_file with regulatory elements.
    For use in multiprocessing.

    Args:
        bed_file_path (str): Path to file containing 
        chrom_name (str): Name of chromosome    
        resolution (int): Resolution of 
        total_parts (int): Total number of parts to divide promoter_dataframe into.
        part (int): After promoter_dataframe_ is divised into X smaller parts, this indicates the part that this process should handle.

    Returns:
        _type_: _description_
    """

    print(f"Calculating predicted promtoer/enhancer interactions on {chrom_name}, res {resolution}. Part {part} of {total_parts}")

    coordinate_np_type = np.uint32

    counter_np_type = np.uint16

    numpy_array_start_size = 10_000

    coordinate_np = np.empty(shape=(numpy_array_start_size,2), dtype = coordinate_np_type)

    counter_np = np.empty(shape=numpy_array_start_size, dtype = counter_np_type)

    current_numpy_index = 0

    promoter_dataframe, enhancer_dataframe = extract_pls_els_from_bed(bed_file_path, split = True)

    promoter_dataframe = promoter_dataframe[promoter_dataframe['chrom'] == chrom_name].reset_index(drop=True)
    enhancer_dataframe = enhancer_dataframe[enhancer_dataframe['chrom'] == chrom_name]

    promoter_dataframe_size = len(promoter_dataframe)

    if total_parts > promoter_dataframe_size:
        raise ValueError(f"Cannot split work into more parts than the number of promoters. Total_parts = {total_parts}, promoter dataframe size: {promoter_dataframe_size}")

    promoters_per_part = round(promoter_dataframe_size / total_parts,0) 


    start_of_part = int(promoters_per_part * part - promoters_per_part)
    # * Make sure we don't overlap with previous part.
    if part != 1: start_of_part += 1

    end_of_part = int(promoters_per_part * part)
    # * if final part, get the rest
    if part == total_parts:
        end_of_part = promoter_dataframe_size - 1


    # * Filter dataframes to only include entries with the chrom we are looking at currently
    promoter_dataframe = promoter_dataframe[start_of_part:end_of_part + 1]


    lowest_prom_start = np.min(np.array(promoter_dataframe['chromStart']))
    highest_prom_start = np.max(np.array(promoter_dataframe['chromStart']))
    lowest_prom_end = np.min(np.array(promoter_dataframe['chromEnd']))
    highest_prom_end = np.max(np.array(promoter_dataframe['chromEnd']))

    enhancer_dataframe = enhancer_dataframe.loc[enhancer_dataframe['chromStart'] <= highest_prom_end + max_distance]
    enhancer_dataframe = enhancer_dataframe.loc[enhancer_dataframe['chromEnd'] >= lowest_prom_start - max_distance]
    # * Make sure dataframes are sorted in the correct order
    promoter_dataframe = promoter_dataframe.sort_values(by=['chromStart','chromEnd'])
    enhancer_dataframe = enhancer_dataframe.sort_values(by=['chromStart','chromEnd'])

    lowest_enh_start = np.min(np.array(enhancer_dataframe['chromStart']))
    highest_enh_start = np.max(np.array(enhancer_dataframe['chromStart']))
    lowest_enh_end = np.min(np.array(enhancer_dataframe['chromEnd']))
    highest_enh_end = np.max(np.array(enhancer_dataframe['chromEnd']))

    # * Iterate through each promoter (only once, so we don't get duplicate counts)
    for prom_row in promoter_dataframe.itertuples():
        prom_start = prom_row[2]; 
        prom_end = prom_row[3]

        # * Round numbers to find exact start and end positions of bins
        # * We find the start of the first bin, and the end of the second bin
        # * Calculate how many bins the promoter spans

        prom_start = int(pelg.round_up_and_down(prom_start, resolution)[0])
        prom_end = int(pelg.round_up_and_down(prom_end, resolution)[1])
        number_of_prom_bins = (prom_end - prom_start) // resolution

        for enh_row in enhancer_dataframe.itertuples():
            enh_start = enh_row[2]; 
            enh_end = enh_row[3]

            # * Round numbers to find exact start and end positions of bins
            enh_start = int(pelg.round_up_and_down(enh_start, resolution)[0])
            enh_end = int(pelg.round_up_and_down(enh_end, resolution)[1])
            
            # * If further apart than 3Mbp, break. This assumes that the enhancer dataframe is sorted by enh_start
            if prom_end < enh_start:
                if abs(enh_start - prom_end) > 3_000_000: 
                    break
            elif enh_end < prom_start:
                if abs(enh_end - prom_start) > 3_000_000: 
                    continue

            # * Calculate how many bins the enhancer spans
            number_of_enh_bins = (enh_end - enh_start) // resolution

            # * For each combination of bins, add a pixel
            for relative_prom_index in range(number_of_prom_bins):
                for relative_enh_index in range(number_of_enh_bins):
                    # * Find current promoter/enhacer bins we are looking at. Should end up divisible with resolution
                    prom_bin = prom_start + relative_prom_index * resolution 
                    enh_bin = enh_start + relative_enh_index * resolution

                    # * Find the indexes of the bins
                    prom_index = prom_bin // resolution 
                    enh_index = enh_bin // resolution 

                    # * Because the matrix is a two dimensional contact matrix, each of the two triangles
                    # * that make up the matrix (if we split along the diagonal) are flipped duplicates of
                    # * each other. Therefore we can choose only one of them to keep count of. 
                    # * If a point appears on the opposite triangle, we instead flip the indexes
                    # * so all points are on the same triangle. 
                    indexes = [prom_index, enh_index]
                    if enh_index < prom_index:
                        indexes = [enh_index, prom_index]

                    #print(chrom_name,total_parts,part,lowest_prom_start,highest_prom_end,lowest_enh_start,highest_enh_start,"\n",prom_start,prom_end,enh_start,enh_end,relative_prom_index,relative_enh_index) #! remove

                    # * Convert to numpy array before comparison to avoid runtimewarning: invalid value encountered in long_scalars
                    indexes = np.array(indexes,dtype=coordinate_np_type)
                    coordinate_index = False
                    for i, e in enumerate(coordinate_np):
                        if np.all(e == indexes):
                            coordinate_index = i
                            break;

                    if coordinate_index:
                        # * If coordinate already exists in array, simply update the counter for the coordinate
                        counter_np[coordinate_index] = counter_np[coordinate_index] + 1
                    else:
                        if current_numpy_index >=  counter_np.size:

                            new_rows = np.empty((numpy_array_start_size), dtype=counter_np_type)
                            counter_np = np.append(counter_np, new_rows)


                            new_rows = np.empty((numpy_array_start_size, coordinate_np.shape[1]), dtype=coordinate_np_type)
                            coordinate_np = np.vstack((coordinate_np, new_rows))
                        # * Add values to numpy arrays

                        

                        coordinate_np[current_numpy_index][0] = indexes[0]
                        coordinate_np[current_numpy_index][1] = indexes[1]
                        counter_np[current_numpy_index] = 1
                        
                        # * Increment index by one 
                        current_numpy_index += 1

    try:
        # * Resize coordinate_np and counter_np to save memory. New size: current_numpy_index 
        coordinate_np = coordinate_np[:current_numpy_index]
        counter_np = counter_np[:current_numpy_index]
        total_bins_with_regulatory_interaction = len(counter_np)


        summed_count = np.sum(counter_np,dtype=np.uint)


    except RuntimeWarning:
        print("Runtimewarning in:", chrom_name, resolution )
    
    dictionary = {'chrom':[chrom_name],
                'resolution':[resolution],
                'total_bins_with_pls_and_els':[total_bins_with_regulatory_interaction],
                'total_count':[summed_count],
                'list_of_counts':[counter_np],
                'list_of_indexes':[coordinate_np]}

    
    # * Convert dict to dataframe
    dataframe = pd.DataFrame.from_dict(dictionary)
    
    print(f"RETURNING {_predict.__name__} with promoter/enhancer bins in {chrom_name}, res {resolution}, part {part} of {total_parts}. Dataframe size: {len(dataframe)}")
    return dataframe

def calculate_promoter_enhancer_bins_multiprocess(bed_file_path : str,
                                    chrom_name : str,
                                    chrom_size : int,
                                    resolution : int,
                                    workers : int,
                                    total_parts : int,
                                    out_file_name : str = False,
                                    ):    
    """_summary_

    Args:
        bed_file_path (str): _description_
        chrom_name (str): _description_
        chrom_size (str): _description_
        resolution (int): _description_
        workers (int): _description_
        total_parts (int): how many parts to divide the task into. Each worker will handle one part at a time.
        out_file_name (str, optional): _description_. Defaults to False.
    """


    print(f"Calculating promoter/enhancer bins in {chrom_name}, res {resolution}")

    dataframe_results : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)
    print(f'Attempting to start subprocesses')

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:

        futures = []
        current_part : int = 1
        while current_part <= total_parts:
            print(f"Submitting one process with function {_predict.__name__} and args {total_parts}, {current_part}.")
            futures.append(executor.submit(_predict, bed_file_path, chrom_name, resolution, total_parts, current_part))
            current_part += 1
        
        print(f'All processes submitted.')
        print(f'Waiting for all processes to finish.')
        print(f'Number of processes: {len(futures)}')

        concurrent.futures.wait(futures)
        print(f'All processes done.')
        
        print(f'Concatenating all dataframes.')
        for f in futures:
            try:
                returned_df = f.result()
                dataframe_results = pd.concat([dataframe_results,returned_df],ignore_index=True)
            except Exception:
                continue
                
        counter_np = np.empty(shape=0, dtype = np.uint16)
        coordinate_np = np.empty(shape=(0,2), dtype = np.uint32)

        for e in np.array(dataframe_results['list_of_counts']):
            counter_np = np.concatenate([counter_np,e])
        for e in np.array(dataframe_results['list_of_indexes']):
            coordinate_np = np.concatenate([coordinate_np,e])
        
        total_bins_in_chrom = chrom_size // resolution
        total_bins_with_regulatory_interaction = np.sum(np.array(dataframe_results['total_bins_with_pls_and_els']))
        min_count = np.min(counter_np)
        max_count = np.max(counter_np)
        average_count = np.mean(counter_np)
        median_count = np.median(counter_np)
        standarddeviation = np.std(counter_np)
        total_count = np.sum(counter_np)

        summed_count = np.sum(np.array(dataframe_results['total_count']))

        dictionary = {'chrom':[chrom_name],
        'resolution':[resolution],
        'total_bins_in_chrom':[total_bins_in_chrom],
        'total_bins_with_pls_and_els':[total_bins_with_regulatory_interaction],
        'min_count':[min_count],
        'max_count':[max_count],
        'average_count':[average_count],
        'median_count':[median_count],
        'standarddeviation':[standarddeviation],
        'total_count':[summed_count],
        'list_of_counts':[counter_np],
        'list_of_indexes':[coordinate_np]}

        return_dataframe = pd.DataFrame(dictionary)

        if not out_file_name:
            out_file_name = f'predicted_chromosome_data_{chrom_name}_{resolution}.csv'

        save_dataframe(dataframe=return_dataframe,
                    file_path=out_file_name,numpy_columns=["list_of_counts","list_of_indexes"])



        print("Dataframe stats:")
        print("Len:",len(return_dataframe))
        print(return_dataframe)

        return dataframe_results



def _predict2(bed_file_path : str, chrom_name : str, resolution : int):

    promoter_dataframe, enhancer_dataframe = extract_pls_els_from_bed(bed_file_path,split=True)

    promoter_dataframe.loc[promoter_dataframe['chrom'] == chrom_name] 
    enhancer_dataframe.loc[enhancer_dataframe['chrom'] == chrom_name] 

    # * Make sure df is sorted
    
    promoter_dataframe = promoter_dataframe.sort_values(by=["chromStart","chromEnd"])
    enhancer_dataframe = enhancer_dataframe.sort_values(by=["chromStart","chromEnd"])



    promoter_indexes = np.zeros(1000)
    current_index = 0

    print("Finding all promoters")

    for row in promoter_dataframe.itertuples():
        prom_start = int(row[2])
        promoter_bin_start = prom_start // resolution
        prom_end = int(row[3])
        promoter_bin_end = prom_end // resolution

        for n in range(promoter_bin_start,promoter_bin_end + 1):
            if current_index == len(promoter_indexes):
                promoter_indexes = np.resize(promoter_indexes, len(promoter_indexes) + 1000)

            promoter_indexes[current_index] = n
            current_index += 1

    enhancer_indexes = np.zeros(1000)
    current_index = 0

    print("Finding all enhancers")

    for row in enhancer_dataframe.itertuples():
        enhancer_start = int(row[2])
        enhancer_bin_start = enhancer_start // resolution
        enhancer_end = int(row[3])
        enhancer_bin_end = enhancer_end // resolution

        for n in range(enhancer_bin_start,enhancer_bin_end + 1):
            if current_index == len(enhancer_indexes):
                enhancer_indexes = np.resize(enhancer_indexes, len(enhancer_indexes) + 1000)
            enhancer_indexes[current_index] = n
            current_index += 1

    
    promoter_indexes = np.resize(promoter_indexes, np.count_nonzero(promoter_indexes))
    enhancer_indexes = np.resize(enhancer_indexes, np.count_nonzero(enhancer_indexes))

    unique_promoter_indexes = np.unique(promoter_indexes)
    unique_enhancer_indexes = np.unique(enhancer_indexes)


    print("Finding all indexes")

    all_indexes = np.empty(shape=(len(promoter_indexes) * len(enhancer_indexes),2),dtype=np.uint16)
    all_indexes_index = 0
    for result in itertools.product(unique_promoter_indexes,unique_enhancer_indexes):
        r1 = int(result[0])
        r2 = int(result[1])
        all_indexes[all_indexes_index] = np.array([r1,r2],dtype=np.uint16)
        all_indexes_index += 1


    unique_elements, counts = np.unique(all_indexes,return_counts=True,axis=1)


    unique_dict = {}
    for element, count in zip(unique_elements, counts):
        unique_dict[element] = count


    print(unique_elements)

    total = np.sum(list(unique_dict.values()))
    bins = len(unique_dict)
    print(total, bins)


    

