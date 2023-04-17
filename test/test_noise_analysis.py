import pandas
import numpy as np
import os
import sys
import time

from line_profiler import LineProfiler #type:ignore

sys.path.append('../src/regelem')
import noise_analysis as na
import files

cooler_path = "../input/H1-hESC.mcool::/resolutions/5000"
output_csv_path = "./output/noise_dataframe_5000.csv"
output_plot_path = "./output/plots/5000_noise_plot.png"
output_plot_path_1000 = "./output/plots/1000_noise_plot.png"

noise_dataframe_path_5000 = "./output/noise_dataframe_h1_5000.csv"
noise_dataframe_path_1000 = "./output/noise_dataframe_h1_1000.csv"
reg_dataframe_path = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"

res = 5000
# h1_dataframe = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
# cooler_file = "../input/H1-hESC.mcool::/resolutions/5000"


def main():
    #na.generate_bias("../input/H1-hESC.mcool::/resolutions/5000",output_path="bias_maybe.tsv")

    # noise_dataframe = na.funct(reg_dataframe_path,cooler_path,res,output_csv_path)

    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)
    lesse_res_dataframe = files.load_dataframe(reg_dataframe_path)
    na.plot_noise(noise_dataframe,lesse_res_dataframe,output_plot_path_1000)
    # #profile_funct()
    # exit()

def profile_funct():
    # Create a new LineProfiler object
    profiler = LineProfiler()

    # Add the function to be profiled to the profiler
    profiler.add_function(na.funct)

    # Run the profiler on the function with arguments
    profiler.run('na.funct(reg_dataframe_path,cooler_path,res)')

    # Print the results
    profiler.print_stats()

if __name__ == "__main__":
    main()