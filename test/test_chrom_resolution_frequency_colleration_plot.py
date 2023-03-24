import os
import pandas as pd
import numpy as np

import sys
sys.path.append('../src/regelem')
import chrom_resolution_frequency_correlation_plot as crfcp

input_folder : str = "../input/predicted/"

def main():
    chrom_and_resolution_average_dataframe = get_output_file_names()

    feature1 = np.array([])
    feature2 = np.array([])
    feature3 = np.array([])
    feature4 = np.array([])

    for row in chrom_and_resolution_average_dataframe.itertuples():
        file_name = input_folder + row[3]
        
        

        file_data = pd.read_csv(file_name,index_col=0)
        chrom_name = file_data['chrom'][0]
        resolution = file_data['resolution'][0]
        average_count = file_data['average_count'][0]
        average_count_per_pls = file_data['total_count'][0] / file_data['total_bins_with_pls_and_els'][0]

        #print(chrom_name, resolution, average_count)

        feature1 = np.append(feature1,chrom_name)
        feature2 = np.append(feature2,resolution)
        feature3 = np.append(feature3,average_count)

        #print(file_data)
        feature4 = np.append(feature4, average_count_per_pls)

    

    crfcp.plot_tripple_correlation(
        feature1, feature2, feature3, 
        "chromosome","resolution","average count",
        file_name="example_tipple_plot.png",
        show = True
    )
    # print(feature1, feature2, feature3)
    # print(feature4)


def get_output_file_names():
    files = os.scandir(input_folder)
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



if __name__ == "__main__":
    main()