import sys
import os
from os.path import dirname, basename, isfile, join

#sys.path.insert(0, '..')
from Constants import *
import numpy as np
import pandas as pd
import PromoterEnhancerListGenerator as pelg
import Plotter as plotter
import matplotlib.pyplot as plt 
import cooler
import time
import math
import psutil
from generate_plots import *
from datetime import datetime
from files import *
import threading
from threading import Thread
import concurrent.futures


test_folder = "./testInput/"
input_folder = "../input/"
output_folder = "../output/"

temp_files = "./temp/"

def main():
    #test_prediction_of_promoter_enhancer()
    # print(OUTPUT_FOLDER)
    #predict_average_promoter_enhancer_count_for_all_chroms_and_resolutions()
    pass
    

def predict_average_promoter_enhancer_count_for_all_chroms_and_resolutions():
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

    print(f"Starting {max_workers} processes with function {calculate_promoter_enhancer_bins.__name__}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for row in chrom_sizes_dataframe.itertuples():
            chrom_name = row[1]
            size = row[2]
            for resolution in resolutions:
                out_file_name = f'./output/predicted_chromosome_data_{chrom_name}_{resolution}.csv'
                if os.path.isfile(out_file_name):
                    print(f'File at path {out_file_name} already exists. Skipping.')
                    continue

                print(f"Submitting one thread with function {calculate_promoter_enhancer_bins.__name__} and args promoter_dataframe, enhancer_dataframe, {chrom_name}, {size}, {resolution}.")
                futures.append(executor.submit(calculate_promoter_enhancer_bins, 
                                                promoter_dataframe, 
                                                enhancer_dataframe,
                                                chrom_name, size, resolution, out_file_name))

        print(f'All processes submitted.')
        print(f'Waiting for all processes to finish.')
        concurrent.futures.wait(futures)

        print(f'All processes done.')
        
        
        columns = futures[0].result().columns
        summary_dataframe = pd.DataFrame(columns=columns)

        for f in futures:
            summary_dataframe = pd.concat([summary_dataframe,f.result()])
        summary_dataframe = summary_dataframe.sort_values(['chrom','resolution']).reindex()

        save_dataframe(summary_dataframe,output_folder+"predicted_chromosome_stats_global.csv")

def _calculate(bed_file_path : str, chrom_name : str, resolution : int, total_parts : int, part : int, max_distance : int = 3_000_000):
    """Calculate possible promoter enhancer interactions based on data from a bed_file with regulatory elements.
    For use in multiprocessing.

    Args:
        bed_file_path (str): _description_
        chrom_name (str): _description_
        resolution (int): _description_
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

    #print("SIZE:", len(promoter_dataframe))


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
            
            #if part == 5: print(part, prom_start, enh_start, )
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
    
    print(f"RETURNING {_calculate.__name__} with promoter/enhancer bins in {chrom_name}, res {resolution}, part {part} of {total_parts}. Dataframe size: {len(dataframe)}")
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
            #bed_file_path = "../input/" + "H1-hESC.7group.bed"
            print(f"Submitting one process with function {_calculate.__name__} and args {total_parts}, {current_part}.")
            futures.append(executor.submit(_calculate, bed_file_path, chrom_name, resolution, total_parts, current_part))
            current_part += 1

        print(f'All processes submitted.')
        print(f'Waiting for all processes to finish.')
        print(f'Number of processes: {len(futures)}')
        concurrent.futures.wait(futures)
        print(f'All processes done.')
        
        print(f'Concatenating all dataframes.')
        for f in futures:
            print(f)
            returned_df = f.result()
            dataframe_results = pd.concat([dataframe_results,returned_df],ignore_index=True)
        
        print(dataframe_results)
        
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

        save_dataframe( dataframe=return_dataframe,
                    filepath=out_file_name)



        print("Dataframe stats:")
        print("Len:",len(return_dataframe))
        print(return_dataframe)

        return dataframe_results
#!
    # * Index counter indicating where to add data in both arrays


def calculate_promoter_enhancer_bins(promoter_dataframe : pd.DataFrame,
                                    enhancer_dataframe : pd.DataFrame,
                                    chrom_name : str,
                                    chrom_size : str,
                                    resolution : int,
                                    out_file_name : str = False,):    

    print(f"Calculating promoter/enhancer bins in {chrom_name}, res {resolution}")

    number_of_rows = number_of_columns = math.ceil(chrom_size / resolution)
    number_of_bins = number_of_rows * number_of_columns

    matrix_shape = (number_of_columns, number_of_rows)

    
    # Filter dataframes to only include entries with the chrom we are looking at currently
    promoter_dataframe = promoter_dataframe[promoter_dataframe['chrom'] == chrom_name]
    enhancer_dataframe = enhancer_dataframe[enhancer_dataframe['chrom'] == chrom_name]

    # Make sure dataframes are sorted in the correct order
    promoter_dataframe.sort_values(by=['chromStart','chromEnd'])
    enhancer_dataframe.sort_values(by=['chromStart','chromEnd'])

    
    if   number_of_rows >= 4_294_967_295:
        print("Using np.uint64 as coordinate datatype. This should never happen, because no chromosome is big enough. What are you doing?")
        coordinate_np_type = np.uint64
    elif number_of_rows >= 65_535:
        print("Using np.uint32 as coordinate datatype.")
        coordinate_np_type = np.uint32
    else:
        print("Using np.uint16 as coordinate datatype.")
        coordinate_np_type = np.uint16

    # 
    counter_np_type = np.uint16
    coordinate_np_type = np.uint16

    numpy_array_start_size = 10_000

    coordinate_np = np.empty(shape=(numpy_array_start_size,2), dtype = coordinate_np_type)

    counter_np = np.empty(shape=numpy_array_start_size, dtype=np.uint)
    # * Index counter indicating where to add data
    current_numpy_index = 0

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
            if abs(enh_start - prom_end) > 3_000_000: 
                break

            # Calculate how many bins the enhoter spans
            number_of_enh_bins = (enh_end - enh_start) // resolution

            # For each combination of bins, add a pixel
            for pb in range(number_of_prom_bins):
                for eb in range(number_of_enh_bins):
                    # Current promoter/enhacer bins we are looking at
                    prom_bin = prom_start + pb * resolution 
                    enh_bin = enh_start + eb * resolution

                    # Find the indexes of the bins
                    prom_index = prom_bin // resolution 
                    enh_index = enh_bin // resolution 

                    # Because the matrix is a two dimensional contact matrix, each of the two triangles
                    # that make up the matrix (if we split along the diagonal) are flipped duplicates of
                    # each other. Therefore we can choose only one of them to keep count of. 
                    # If a point appears on the opposite triangle, we instead flip the indexes
                    # so all points are on the same triangle. 
                    indexes = [prom_index, enh_index]
                    if enh_index < prom_index:
                        indexes = [enh_index, prom_index]

                    # Convert to numpy array before comparison to avoid runtimewarning: invalid value encountered in long_scalars
                    indexes = np.array(indexes,dtype=coordinate_np_type)
                    coordinate_index = False
                    for i, e in enumerate(coordinate_np):
                        if np.all(e == indexes):
                            coordinate_index = i
                            break;

                    if coordinate_index:
                        # If coordinate already exists in array, simply update the counter for the coordinate
                        counter_np[coordinate_index] = counter_np[coordinate_index] + 1
                    else:
                        if current_numpy_index >=  counter_np.size:

                            new_rows = np.empty((numpy_array_start_size), dtype=counter_np_type)
                            counter_np = np.append(counter_np, new_rows)


                            new_rows = np.empty((numpy_array_start_size, coordinate_np.shape[1]), dtype=coordinate_np_type)
                            coordinate_np = np.vstack((coordinate_np, new_rows))
                        # Add values to numpy arrays

                        coordinate_np[current_numpy_index][0] = indexes[0]
                        coordinate_np[current_numpy_index][1] = indexes[1]
                        counter_np[current_numpy_index] = 1
                        
                        # Increment index by one 
                        current_numpy_index += 1

                        
                        #coordinate_np = np.append(coordinate_np,np.array([indexes]),axis=0)
                        #counter_np = np.append(counter_np,1)

    try:
        # Resize coordinate_np and counter_np to save memory. New size: current_numpy_index 
        coordinate_np = coordinate_np[:current_numpy_index]
        counter_np = counter_np[:current_numpy_index]
        total_bins_with_regulatory_interaction = len(counter_np)


        summed_count = np.sum(counter_np,dtype=np.uint)
        min_count = np.min(counter_np)
        max_count = np.min(counter_np)
        median_count = np.median(counter_np)
        std_count = np.std(counter_np)
        average_non_zero = summed_count / total_bins_with_regulatory_interaction
        average_regulatory_count_per_bin_total = summed_count / number_of_bins

    except RuntimeWarning:
        print("Runtimewarning in:", chrom_name, resolution )
    
    dictionary = {'chrom':[chrom_name],'resolution':[resolution],'total_bins_in_chrom':[number_of_bins],
                'total_bins_with_pls_and_els':[total_bins_with_regulatory_interaction],
                'min_count':[min_count],'max_count':[max_count],'average_count':[average_regulatory_count_per_bin_total],
                'median_count':[median_count],'standarddeviation':[std_count],'total_count':[summed_count],
                'list_of_counts':[counter_np],'list_of_indexes':[coordinate_np]}
    
    # Convert dict to dataframe
    dataframe = pd.DataFrame.from_dict(dictionary)

    if not out_file_name:
        out_file_name = f'./output/predicted_chromosome_data_{chrom_name}_{resolution}.csv'

    save_dataframe( dataframe=dataframe,
                    filepath=out_file_name)
    
    print(f"RETURNING promoter/enhancer bins in {chrom_name}, res {resolution}")
    return dataframe

if __name__ == "__main__":
    main()

