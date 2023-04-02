import warnings
#from generate_plots import *
import generate_plots
import os
import pandas as pd
import numpy as np
import Plotter as plotter
import matplotlib.pyplot as plt #type:ignore

output_folder = "../../input/predicted/"


chrom_sizes = {
"chr1":    248956422,
"chr2":    242193529,
"chr3":    198295559,
"chr4":    190214555,
"chr5":    181538259,
"chr6":    170805979,
"chr7":    159345973,
"chrX":    156040895,
"chr8":    145138636,
"chr9":    138394717,
"chr11":   135086622,
"chr10":   133797422,
"chr12":   133275309,
"chr13":   114364328,
"chr14":   107043718,
"chr15":   101991189,
"chr16":   90338345,
"chr17":   83257441,
"chr18":   80373285,
"chr20":   64444167,
"chr19":   58617616,
"chrY":   57227415,
"chr22":   50818468,
"chr21":   46709983}

def bins_within_diagonal(region_len : int, resolution : int, diag_len : int):
    bins_on_diagonal = region_len // resolution
    diag_len_bins = diag_len // resolution
    one_diag = (diag_len * 2 + 1) // resolution

    bins_before_diag_len = (region_len - diag_len) // resolution

    bins_within_diag = bins_before_diag_len * one_diag

    for i in range(diag_len_bins + 1):
        bins_on_diagonal2 = diag_len_bins - i
        one_diag = (bins_on_diagonal2 * 2 + 1)
        bins_within_diag += one_diag

    return bins_within_diag


def main():
    plot_average_regulatory_interaction()

def get_output_files_and_collect_to_dataframe():
    files = os.scandir(output_folder)
    dictionary = {"chrom_name":[],"resolution":[],"file_name":[]}

    for entry in files:
        file_name = entry.name
        if "predicted_chromosome_data_" == file_name[:26]:
            file_name_split = file_name.split("_")
            chrom_name = file_name_split[3]
            resolution = file_name_split[4].split(".")[0]

            dictionary["chrom_name"].append(chrom_name)
            dictionary["resolution"].append(resolution)
            dictionary["file_name"].append(file_name)

    dataframe = pd.DataFrame(dictionary)
    return dataframe

def plot_average_regulatory_interaction():
    def _read_file_data(file_name : str):
        file_path = output_folder + file_name

        data = pd.read_csv(file_path)

        data = data.drop('Unnamed: 0',axis=1)
        return data
        
    #* Collect all the data from each chromosome and resolution into a single dataframe
    dataframe = get_output_files_and_collect_to_dataframe()

    # * Get all chromosome names
    chrom_names = pd.unique(dataframe["chrom_name"])

    # * Find out which resolutions are available for all chromosomes. If any are missing, we don't bother plotting them.
    resolutions_in_all_chroms = np.array([])
    for name in chrom_names:
        if resolutions_in_all_chroms.size == 0:
            resolutions_in_all_chroms = np.array(dataframe[dataframe["chrom_name"] == name]["resolution"])
        else:
            available_resolutions = np.array(dataframe[dataframe["chrom_name"] == name]["resolution"])

            for res in resolutions_in_all_chroms:

                if res not in available_resolutions:
                    resolutions_in_all_chroms = resolutions_in_all_chroms[resolutions_in_all_chroms != res]

    # * Remove any rows that have resoltuions that arent available for all chromosomes.
    reduced_dataframe = pd.DataFrame(columns=dataframe.columns)
    for res in resolutions_in_all_chroms:
        part_df = dataframe[dataframe["resolution"] == res]
        reduced_dataframe = pd.concat([reduced_dataframe,part_df])
    
    # * Dictionary to turn into dataframe later
    chr_wide_dict = {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}

    # * Collect data from all chromsomes for each resolution and make a dataframe
    for res in resolutions_in_all_chroms:
        average_regulatory_count_per_bin_total = 0
        average_regulatory_count_per_pe_bin_total = 0
        total_max_count = 0
        highest_max_count = 0

        chr_wide_dict["resolution"].append(int(res))
        for name in chrom_names:
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("ignore", category=UserWarning)
                row = reduced_dataframe[reduced_dataframe['chrom_name'] == name][reduced_dataframe['resolution'] == res]
                rowArray = np.array(row)[0]
                
                file_name = rowArray[2]
                data = _read_file_data(file_name)
                
                average_regulatory_count_per_bin_total += int(data['average_count'])
                average_regulatory_count_per_pe_bin_total += int(data['total_count']) / int(data['total_bins_with_pls_and_els'])
                
                total_max_count += int(data['max_count'])
                if int(data['max_count']) > highest_max_count:
                    highest_max_count = int(data['max_count'])

        average_max_count = total_max_count / len(chrom_names)
        average_regulatory_count_per_bin = average_regulatory_count_per_bin_total / len(chrom_names)
        average_regulatory_count_per_pe_bin = average_regulatory_count_per_pe_bin_total / len(chrom_names)

        chr_wide_dict["max_count"].append(highest_max_count)
        chr_wide_dict["average_regulatory_count_per_bin"].append(average_regulatory_count_per_bin)
        chr_wide_dict["average_regulatory_count_per_pe_bin"].append(average_regulatory_count_per_pe_bin)
    
    
    chr_wide_dataframe = pd.DataFrame(chr_wide_dict).sort_values(['resolution']).reset_index(drop=True)
    chr_wide_dataframe.drop(chr_wide_dataframe.tail(1).index,inplace=True)

    coolerPlotter = plotter.CoolerPlotter()

    resolution = np.array(chr_wide_dataframe['resolution'])
    average_count_per_pe_bin = chr_wide_dataframe['average_regulatory_count_per_pe_bin']
    average_count_per_bin = chr_wide_dataframe['average_regulatory_count_per_bin']

    # * Get data from files
    modle_dictionary_5000 = reformat_statistics_file_and_get_dict("../../input/dataframes/modle_binsize5000/statistics/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.csv", 5000)
    modle_dictionary_1000 = reformat_statistics_file_and_get_dict("../../input/dataframes/modle_binsize1000/statistics/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv",1000)


    H1_dictionary_5000 = reformat_statistics_file_and_get_dict("../../input/dataframes/H1_binsize5000/statistics/5000_H1-hESC.7group.bed", 5000)
    H1_dictionary_1000 = reformat_statistics_file_and_get_dict("../../input/dataframes/H1_binsize1000/statistics/1000_H1-hESC.7group.bed", 1000)

    # * Calculate
    modle_average_regulatory_count_per_pe_bin = np.concatenate((modle_dictionary_1000['average_regulatory_count_per_pe_bin'],modle_dictionary_5000['average_regulatory_count_per_pe_bin']))
    modle_average_regulatory_count_per_bin = np.concatenate((modle_dictionary_1000['average_regulatory_count_per_bin'],modle_dictionary_5000['average_regulatory_count_per_bin']))

    modle_res = np.concatenate((modle_dictionary_1000['resolution'],modle_dictionary_5000['resolution']))
    modle_max_count = np.concatenate((modle_dictionary_1000['max_count'],modle_dictionary_5000['max_count']))
    
    H1_res = np.concatenate((H1_dictionary_1000['resolution'],H1_dictionary_5000['resolution']))
    H1_max_count = np.concatenate((H1_dictionary_1000['max_count'],H1_dictionary_5000['max_count']))
    H1_average_regulatory_count_per_pe_bin = np.concatenate((H1_dictionary_1000['average_regulatory_count_per_pe_bin'],H1_dictionary_5000['average_regulatory_count_per_pe_bin']))
    H1_average_regulator
    y_count_per_bin = np.concatenate((H1_dictionary_1000['average_regulatory_count_per_bin'],H1_dictionary_5000['average_regulatory_count_per_bin']))


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(40, 15))

    # * Average per PE bin: Total counts divided by number of bins with at least one promoter and enhancer pair
    ax1.plot(resolution.astype('str'), average_count_per_pe_bin, marker='o',
            label="Average number of possible regulatory interaction per bin with at least 1 promoter and 1 enhancer, genomewide.")
    ax1.plot(modle_res, modle_average_regulatory_count_per_pe_bin, marker='o',
            label="Average number of regulatory interaction simulated by modle per bin with at least 1 promoter and 1 enhancer, genomewide.")
    ax1.plot(H1_res, H1_average_regulatory_count_per_pe_bin, marker='o',
            label="Average number of regulatory interaction registered in H1-hESC per bin with at least 1 promoter and 1 enhancer, genomewide.")

    ax1.set_xlabel("Resolution")
    ax1.set_ylabel("Average contact per bin")
    ax1.set_yticks([0, 1, 5, 20, 40])
    ax1.set_xticks(resolution.astype('str'))
    ax1.set_yscale('log')
    ax1.grid()
    ax1.legend()

    # * Average per bin: Total counts divided by number of bins in chromosome. Calculation done on each chromosome, then collected and divided by number of chromosomes. 
    ax2.plot(resolution.astype('str'), average_count_per_bin, marker='o',
            label="Average number of possible regulatory interaction per bin, genomewide.")
    ax2.plot(modle_res, modle_average_regulatory_count_per_bin, marker='o',
            label="Average number of regulatory interaction per bin simulated by modle, genomewide.")
    ax2.plot(H1_res, H1_average_regulatory_count_per_bin, marker='o',
            label="Average number of regulatory interaction per bin registered in H1-hESC, genomewide.")

    ax2.set_xlabel("Resolution")
    ax2.set_ylabel("Average contact per bin")
    ax2.grid()
    ax2.legend()

    # * Save and show figure
    fig.show()
    fig.savefig("combined_plots.png")
    coolerPlotter.view_plot("combined_plots.png")

def reformat_statistics_file_and_get_dict(file_path : str, resolution : int) -> dict:
    """Reads a dataframe from a file and extracts the 'max_count', 'total_count' and 'total_bins_with_pls_and_els' columns. 


    Args:
        file_path (str): path to dataframe file (.csv)
        resolution (int): resolution value. To be added to dictionary

    Returns:
        dict: dictionary {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}
    """
    resolution = resolution

    dataframe = pd.read_csv(file_path)

    dictionary = {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}


    max_count = np.mean(dataframe['max_count'])
    average_regulatory_count_per_pe_bin = np.sum(dataframe['total_count']) / np.sum(dataframe['total_bins_with_pls_and_els'])


    # * Find all bins without the diagonal
    total_bins = np.sum([bins_within_diagonal(size,resolution,3_000_000) for size in chrom_sizes.values()])
    average_regulatory_count_per_bin = np.sum(dataframe['total_count']) / total_bins

    dictionary['max_count'].append(max_count) 
    dictionary['average_regulatory_count_per_bin'].append(average_regulatory_count_per_bin)
    dictionary['average_regulatory_count_per_pe_bin'].append(average_regulatory_count_per_pe_bin) 
    dictionary['resolution'].append(str(resolution))
    
    return dictionary
    


if __name__ == "__main__":
    main()