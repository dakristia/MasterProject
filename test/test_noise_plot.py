import pandas as pd
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt #type:ignore

import noise_analysis as na
import files
import plot_utils

def plot_noise(noise_dataframe : pd.DataFrame, lesser_res_dataframe : pd.DataFrame, output_path = False):
    ""

    #TODO: Move this function to test_noise_analysis.py

    noise_dataframe = pd.DataFrame(noise_dataframe).sort_values(by=["chrom","bin1_start","bin2_start"])

    # Remove duplicate pixels. This will remove some pairs, but we only care about specific pixels right now. 
    duplicates = noise_dataframe.loc[:,["chrom","bin1_start","bin1_end","bin2_start","bin2_end"]].duplicated(keep="first",)
    
    noise_dataframe = noise_dataframe[~duplicates].sort_values(by=["chrom","bin1_start","bin2_start"])

    chrom_array = np.array(noise_dataframe["chrom"])
    count_array = np.array(noise_dataframe["count"])
    noise_array = np.array(noise_dataframe["noise"])
    noise_count_ratio_array = noise_array / count_array

    adjacent_columns = ["down","downleft","left","upleft","up","upright","right","downright"]

    
    max_noise_array = np.array(noise_dataframe[adjacent_columns].max(axis=1))
    min_noise_array = np.array(noise_dataframe[adjacent_columns].min(axis=1))
    median_noise_array = np.array(noise_dataframe[adjacent_columns].median(axis=1))
    average_noise_array = np.array(noise_dataframe[adjacent_columns].mean(axis=1))

    # create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12,10))

    # create a set of unique labels
    unique_labels = set(chrom_array)
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap #type:ignore
    viridis = plt.cm.get_cmap('gist_rainbow', len(unique_labels))
    color_choices = viridis(np.linspace(0, 1, 24))
    color_dict = {label: color_choices[i] for i, label in enumerate(unique_labels)}
    colors = [color_dict[label] for label in chrom_array]
    
    # * Plot 1

    for label in unique_labels:
        indices = [i for i, x in enumerate(chrom_array) if x == label]
        uq_count_array = [count_array[i] for i in indices]
        uq_noise_array = [noise_array[i]/8 for i in indices]
        scatter = axs[0, 0].scatter(uq_count_array, uq_noise_array, alpha = 0.4, label=label ,color=color_dict[label], s = 20)

    axs[0, 0].plot(count_array, count_array, color="black")
    #axs[0, 0].plot(count_array, count_array * 8, color="brown")
    #scatter = axs[0, 0].scatter(count_array, noise_array, c=colors, s = 20)
    axs[0, 0].set_title('Count vs noise per pair')
    axs[0, 0].set_xlabel('Center frequency "count"')
    axs[0, 0].set_ylabel('Average surrounding\nfrequency "noise"')
    #axs[0, 0].set_yscale('log')
    plt.rcParams['legend.fontsize'] = 'small'
    axs[0, 0].legend(loc='upper center', columnspacing=0.2, bbox_to_anchor=(0.22, 1), ncol=3, prop={'size': 6})


    # * Plot 2
    box = axs[0, 1].boxplot([count_array, min_noise_array, max_noise_array, median_noise_array, average_noise_array], labels=["Count","Min\nnoise","Max\nnoise","Median\nnoise","Average\nnoise"])
    axs[0, 1].set_yscale('log')



    # * Plot 3
    label_int = np.unique(chrom_array, return_inverse=True)[1]
    count_sums = np.bincount(label_int, weights=count_array)
    noise_sums = np.bincount(label_int, weights=noise_array)

    median_noises = np.array([])
    noise_lengths = np.array([])

    for label in np.unique(chrom_array):
        label_values = noise_array[chrom_array == label]  # Select the values connected to the current label
        median = np.median(np.sort(label_values))
        median_noises = np.append(median_noises,median)
        noise_lengths = np.append(noise_lengths,len(label_values))

    noise_averages = noise_sums/noise_lengths

    width = 0.25
    x = np.arange(len(count_sums))

    rects1 = axs[1, 0].bar(x - width, count_sums, width, label='Count')
    rects2 = axs[1, 0].bar(x, noise_sums, width, label='Noise')
    #rects3 = axs[1, 0].bar(x + width, noise_averages, width, label='Mean Noise', color = 'purple')

    # add labels for the x and y axes
    axs[1, 0].set_xticks(x)
    axs[1, 0].set_xticklabels(np.unique(chrom_array),rotation=45, ha='right')

    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlabel('Label')
    axs[1, 0].set_ylabel('Sum')

    plt.rcParams['legend.fontsize'] = 'small'

    axs[1, 0].legend(prop={'size': 6})
    # * Plot 4
    ratio_sums = np.bincount(label_int, weights=noise_count_ratio_array)
    count = np.bincount(label_int)
    ratio_average = ratio_sums / count

    width = 0.5

    rects3 = axs[1, 1].bar(x, ratio_average, width, label='Noise per count', color="green")
    
    axs[1, 1].set_xticks(x)
    axs[1, 1].set_xticklabels(np.unique(chrom_array),rotation=45, ha='right')
    axs[1, 1].set_yscale('log')
    axs[1, 1].legend()

    # * Show the plot
    plt.show()
    if output_path: 
        plt.savefig(output_path)
        plot_utils.show_plot(output_path)