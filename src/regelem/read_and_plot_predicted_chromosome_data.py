import warnings
#from generate_plots import *
import generate_plots
import os
import pandas as pd
import numpy as np
import Plotter as plotter
import matplotlib.pyplot as plt #type:ignore

output_folder = "../../output/predicted/"
def main():
    plot_average_regulatory_interaction()

def get_output_file_names():
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
        
    dataframe = get_output_file_names()

    chrom_names = pd.unique(dataframe["chrom_name"])

    resolutions_in_all_chroms = np.array([])
    for name in chrom_names:
        if resolutions_in_all_chroms.size == 0:
            resolutions_in_all_chroms = np.array(dataframe[dataframe["chrom_name"] == name]["resolution"])
        else:
            available_resolutions = np.array(dataframe[dataframe["chrom_name"] == name]["resolution"])

            for res in resolutions_in_all_chroms:

                if res not in available_resolutions:
                    resolutions_in_all_chroms = resolutions_in_all_chroms[resolutions_in_all_chroms != res]

    reduced_dataframe = pd.DataFrame(columns=dataframe.columns)
    for res in resolutions_in_all_chroms:
        part_df = dataframe[dataframe["resolution"] == res]
        reduced_dataframe = pd.concat([reduced_dataframe,part_df])
    
    chr_wide_dict = {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}

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
    print(chr_wide_dataframe)

    coolerPlotter = plotter.CoolerPlotter()

    resolution = np.array(chr_wide_dataframe['resolution'])
    average_count_per_pe_bin = chr_wide_dataframe['average_regulatory_count_per_pe_bin']

    # #TODO: Make this entire file less tacky
    # modle_dictionary_5000 = reformat_statistics_file_and_get_dict("../../output/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.csv", 5000)
    # modle_dictionary_1000 = reformat_statistics_file_and_get_dict("../../output/old/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv",1000)

    # modle_res = np.concatenate((modle_dictionary_1000['resolution'],modle_dictionary_5000['resolution']))

    # modle_max_count = np.concatenate((modle_dictionary_1000['max_count'],modle_dictionary_5000['max_count']))

    # modle_average_regulatory_count_per_pe_bin = np.concatenate((modle_dictionary_1000['average_regulatory_count_per_pe_bin'],modle_dictionary_5000['average_regulatory_count_per_pe_bin']))

    # H1_dictionary_5000 = reformat_statistics_file_and_get_dict("../../output/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.csv", 5000)
    # H1_dictionary_1000 = reformat_statistics_file_and_get_dict("../../output/old/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv",1000)

    # H1_res = np.concatenate((H1_dictionary_1000['resolution'],H1_dictionary_5000['resolution']))

    # H1_max_count = np.concatenate((H1_dictionary_1000['max_count'],H1_dictionary_5000['max_count']))

    # H1_average_regulatory_count_per_pe_bin = np.concatenate((H1_dictionary_1000['average_regulatory_count_per_pe_bin'],H1_dictionary_5000['average_regulatory_count_per_pe_bin']))

    fig, ax = plt.subplots(figsize=(20,15))

    ax.plot(resolution.astype('str'),average_count_per_pe_bin, 
            label = "Average number of possible regulatory interaction per bin with at least 1 promoter and 1 enhancer, genomewide.")
    # ax.plot(modle_res,modle_average_regulatory_count_per_pe_bin, 
    #         label = "Average number of regulatory interaction simulated by modle per bin with at least 1 promoter and 1 enhancer, genomewide.")
    # ax.plot(H1_res,H1_average_regulatory_count_per_pe_bin, 
    #         label = "Average number of regulatory interaction registered in H1-hESC per bin with at least 1 promoter and 1 enhancer, genomewide.")


    plt.xlabel("Resolution")
    plt.ylabel("Bins")

    # Ensure 1 is present on the y axis
    plt.yticks(list(plt.yticks()[0]) + [1])
    ax.yaxis.set_ticks([0,1,5,20,40])

    #Ensures x ticks are spaced out evenly
    ax.xaxis.set_ticks(resolution.astype('str'))
    #ax.ticklabel_format(useOffset=True, style='plain')
    
    plt.grid()
    plt.legend()

    fig.show()
    fig.savefig("test_line_plot.png")
    coolerPlotter.view_plot("test_line_plot.png")

def reformat_statistics_file_and_get_dict(file_path : str, resolution : int) -> dict:
    """Reads a dataframe from a file and extracts the 'max_count', 'total_count' and 'total_bins_with_pls_and_els' columns. 


    Args:
        file_path (str): path to dataframe file (.csv)
        resolution (int): resolution value. To be added to dictionary

    Returns:
        dict: dictionary {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}
    """
    resolution = str(resolution)

    dataframe = pd.read_csv(file_path)

    dictionary = {"resolution":[],"max_count":[],"average_regulatory_count_per_bin":[], "average_regulatory_count_per_pe_bin":[]}

    max_count = np.mean(dataframe['max_count'])
    average_regulatory_count_per_pe_bin = np.sum(dataframe['total_count']) / np.sum(dataframe['total_bins_with_pls_and_els'])

    dictionary['max_count'].append(max_count) 
    dictionary['average_regulatory_count_per_pe_bin'].append(average_regulatory_count_per_pe_bin) 
    dictionary['resolution'].append(resolution)
    
    return dictionary
    

def add_statistics_file_line(file_paths : np.array, resolutions : np.array, ax):
    print(type(ax))

    for file_path in file_paths:

        dataframe = pd.read_csv(file_path)

        #value = 

if __name__ == "__main__":
    main()