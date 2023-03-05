import sys
import os
from os.path import dirname, basename, isfile, join

sys.path.insert(0, '..')
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
import threading
from threading import Thread
import concurrent.futures


test_folder = "./testInput/"
input_folder = "../input/"
output_folder = "./output/"
temp_files = "./temp/"

def main():
    #test_prediction_of_promoter_enhancer()
    predict_average_promoter_enhancer_count_for_all_chroms_and_resolutions()


def predict_average_promoter_enhancer_count_for_all_chroms_and_resolutions():
    max_workers = 20
    chrom_sizes_file_path = input_folder + "hg38.chrom.sizes"
    bed_file_path = input_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"
    cooler_file_path = input_folder + "4DNFI9GMP2J8.mcool"

    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])

    # Get promoters and enhancers from bed file
    promoter_enhancer_dataframe = pd.read_csv(bed_file_path,delim_whitespace=True,header=None,names=DATAFRAME_COLUMNS_BED)
    # Filter 'type' column to only include PLS, pELS or dELS (makes it easier to check this column later)
    promoter_enhancer_dataframe = pelg.filter_type_in_dataframe(promoter_enhancer_dataframe)
    # Split dataframe into one for promoters and one for enhancers
    promoter_dataframe, enhancer_dataframe = pelg.split_df_to_pls_els(promoter_enhancer_dataframe)

    # Array with all resolutions we want to calculate for
    resolutions = np.array([50,100,500,1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000])

    print(chrom_sizes_dataframe)
    print(f"Starting {max_workers} threads with function {calculate_promoter_enhancer_bins.__name__}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        #futures = [executor.submit(calculate_promoter_enhancer_bins, promoter_dataframe,enhancer_dataframe, res) for res in resolutions]

        futures = []
        for row in chrom_sizes_dataframe.itertuples():
            
            
            for res in resolutions:
                
                chrom = row[1]
                size = row[2]

                print(f"Submitting one thread with function {calculate_promoter_enhancer_bins.__name__} and args promoter_dataframe, enhancer_dataframe, {chrom}, {size}, {res}.")
                futures.append(executor.submit(calculate_promoter_enhancer_bins, 
                                                promoter_dataframe, 
                                                enhancer_dataframe,
                                                chrom, size, res))

        print(f'All threads submitted.')
        print(f'Waiting for all threads to finish.')
        concurrent.futures.wait(futures)

        print(f'All threads done.')
        
        
        columns = futures[0].result().columns
        summary_dataframe = pd.DataFrame(columns=columns)

        for f in futures:
            summary_dataframe = pd.concat([summary_dataframe,f.result()])
        print(summary_dataframe)
        summary_dataframe = summary_dataframe.sort_values(['chrom','resolution']).reindex()

        save_dataframe(summary_dataframe,output_folder+"predicted_chromosome_stats_global.csv")

        



def calculate_promoter_enhancer_bins(promoter_dataframe : pd.DataFrame,
                                    enhancer_dataframe : pd.DataFrame,
                                    chrom_name : str,
                                    chrom_size : str,
                                    resolution : int):    

    number_of_rows = number_of_columns = math.ceil(chrom_size / resolution)
    number_of_bins = number_of_rows * number_of_columns
    matrix_shape = (number_of_columns, number_of_rows)

    
    # Filter dataframes to only include entries with the chrom we are looking at currently
    promoter_dataframe = promoter_dataframe[promoter_dataframe['chrom'] == chrom_name]
    enhancer_dataframe = enhancer_dataframe[enhancer_dataframe['chrom'] == chrom_name]

    # Make sure dataframes are sorted in the correct order
    promoter_dataframe.sort_values(by=['chromStart','chromEnd'])
    enhancer_dataframe.sort_values(by=['chromStart','chromEnd'])

    print(promoter_dataframe)
    print(enhancer_dataframe)

    # Numpy array to contain coordinates we have currently found
    coordinate_np = np.empty(shape=(0,2), dtype=int)
    counter_np = np.empty(shape=0,dtype=int)

    cc = 0

    for prom_row in promoter_dataframe.itertuples():
        prom_start = prom_row[2]; 
        prom_end = prom_row[3]
        
        # Round numbers to find exact start and end positions of bins
        # We find the start of the first bin, and the end of the second bin
        # Calculate how many bins the promoter spans

        prom_start = int(pelg.round_up_and_down(prom_start, resolution)[0])
        prom_end = int(pelg.round_up_and_down(prom_end, resolution)[1])
        number_of_prom_bins = (prom_end - prom_start) // resolution


        for enh_row in enhancer_dataframe.itertuples():
            enh_start = enh_row[2]; 
            enh_end = enh_row[3]

            # Round numbers to find exact start and end positions of bins
            enh_start = int(pelg.round_up_and_down(enh_start, resolution)[0])
            enh_end = int(pelg.round_up_and_down(enh_end, resolution)[1])
            

            # If further apart than 3Mbp, break. This assumes that the enhancer dataframe is sorted by enh_start
            if abs(enh_start - prom_end) > 3_000_000: 
                #print(f"Promoter in {chrom_name} at {prom_start}-{prom_end} and enhancer at {enh_start}-{enh_end} greater than 3Mbp apart. Breaking")
                break

            # Calculate how many bins the enhoter spans
            number_of_enh_bins = (enh_end - enh_start) // resolution

            print(prom_start,prom_end,enh_start,enh_end)
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

                    coordiante_index = False
                    for i, e in enumerate(coordinate_np):
                        if np.all(e == indexes):
                            coordiante_index = i
                            break;
                    # Find out if list already is in coordinate_np. 
                    #already_in = np.all(np.isin(coordinate_np, indexes).reshape(coordinate_np.shape), axis=1)
                    # already_in = np.isin(coordinate_np,[indexes])
                    # print(f"{[indexes]} already in?")
                    # print(already_in)
                    # print(coordinate_np)
                    # If already in coordinate_np, add counter to correlating index in counter_np
                    if coordiante_index:
                        
                        # Create mask for use in np.where
                        #mask = already_in
                        # Find the index of the array that is already in coordinate_np
                        #idx = np.where(mask)[0]
                        #index = idx[0]
                        #prev_value = 
                        counter_np[coordiante_index] = counter_np[coordiante_index] + 1
                        print(f"Value at index {indexes} updated from {counter_np[coordiante_index] - 1} to {counter_np[coordiante_index]}")


                        #cc += 1
                        # if cc == 10:
                        #     exit()
                    else:
                        print("Adding new index:", indexes)
                        coordinate_np = np.append(coordinate_np,np.array([indexes]),axis=0)
                        #print(coordinate_np)
                        counter_np = np.append(counter_np,1)


            # print(f"Length: {len(coordinate_np)}")
            # print(f"Size of coordinates: {sys.getsizeof(coordinate_np) // 1_048_576}MB") 
            # print(f"Size of counter:{sys.getsizeof(counter_np) // 1_048_576}MB")

    total_bins_with_pls_and_els = len(counter_np)
    average_pls_els_per_bin = np.sum(counter_np) / number_of_bins
    average_non_zero = np.sum(counter_np) / total_bins_with_pls_and_els
    number_of_bins_with_promoter_enhancer_bin = len(coordinate_np)
    
    #! TODO: SAVE TO DATAFRAME WITH COLUMNS 
    #! [CHROMOSOME, RESOLUTION, TOTAL_BINS_IN_CHROM, TOTAL_BINS_WITH_PROMOTER_ENHANCER, MIN_COUNT, MAX_COUNT, AVERAGE_COUNT, MEDIAN_COUNT
    #! STANDARDDEVIATION, TOTALCOUNT, LIST_OF_COUNTS, LIST_OF_INDEXES]
    #append_text_to_file(f'Number of bins containing promoter and enhancer for {chrom_name} at {resolution}bp resolution: {len(coordinate_np)}\n{chrom_name} {resolution}',
    #                    f'./output/average_promoter_enhancer_interaction_per_bin_{chrom_name}_{resolution}.txt')
    print(coordinate_np,counter_np)
    # Make dictionary with collected data
    dictionary = {'chrom':[chrom_name],'resolution':[resolution],'total_bins_in_chrom':[number_of_bins],
                'total_bins_with_pls_and_els':[total_bins_with_pls_and_els],
                'min_count':[np.min(counter_np)],'max_count':[np.max(counter_np)],'average_count':[average_pls_els_per_bin],
                'median_count':[np.median(counter_np)],'standarddeviation':[np.std(counter_np)],'total_count':[np.sum(counter_np)],
                'list_of_counts':[counter_np],'list_of_indexes':[coordinate_np]}
    
    # Convert dict to dataframe
    dataframe = pd.DataFrame.from_dict(dictionary)

    save_dataframe( dataframe=dataframe,
                    filename=f'./output/predicted_chromosome_data_{chrom_name}_{resolution}.csv')
    
    return dataframe

if __name__ == "__main__":
    main()

