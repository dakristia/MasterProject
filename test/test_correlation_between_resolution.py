import numpy as np

import sys
sys.path.append('../src/regelem')
import correlation_between_resolutions
import files

def main():
    print("ayo")
#"output/"
    input_folder = "../input/"
    output_folder = "../output/"

    #* Modle FIles
    dataframe1_file_name = "cache/chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.chr1.csv"
    dataframe2_file_name = "chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv"

    #* H1 files genome wide

    dataframe1_file_name = f"{input_folder}dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
    dataframe2_file_name = f"{input_folder}dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    
    # # * Test Files
    # dataframe1_file_name = f"{input_folder}dataframes/H1_binsize5000/5000_H1-hESC.7group.bed.chrY.csv"
    # dataframe2_file_name = f"{input_folder}dataframes/H1_binsize1000/1000_H1-hESC.7group.bed.chrY.csv"

    dataframe1_filepath = dataframe1_file_name
    dataframe2_filepath = dataframe2_file_name

    dataframe1 = files.load_dataframe(dataframe1_filepath)
    dataframe2 = files.load_dataframe(dataframe2_filepath)
    
    dataframe_array = [dataframe1,dataframe2]

    print(dataframe2)

    label1 = "H1 5000bp"
    label2 = "H1 1000bp"
    labels_array = np.array([label1,label2])

    #correlation_between_resolutions.box_plots(dataframe_array,labels_array)

    #correlation_between_resolutions.violin_plot(dataframe_array,labels_array)

    correlation_between_resolutions.plot_correlation_between_resolution(dataframe1,dataframe2,label1=label1,label2=label2,output_path="../../output/plots/correlation_H1_5000_1000.png",res1= 5000,res2=1000,logScale= True ,scatter_plot = True)

if __name__ == "__main__":
    main()