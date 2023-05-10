import pandas as pd
import numpy as np
import cooler
from line_profiler import LineProfiler # type:ignore


import sys
sys.path.append('../src/regelem')
import plot_utils
import files
import Constants

def main():
    remove_duplicate_noise_from_noise_dataframe()


def remove_duplicate_noise_from_noise_dataframe():
    resolution = 1000
    noise_dataframe_path_1000 = "./output/noise_dataframe_h1_1000.csv"
    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)

    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"
    cooler_matrix = cooler.Cooler(cooler_file).matrix(balance=False,as_pixels=True,join=True)

    close_diag_df = noise_dataframe[abs(noise_dataframe['bin1_start'] - noise_dataframe['bin2_start']) == 1000]

    count_and_noise_columns = ["count","down","downleft","left","upleft","up","upright","right","downright"]

    upper_triangle = close_diag_df[close_diag_df['bin2_start'] > close_diag_df['bin1_start']]

    lower_triangle = close_diag_df[close_diag_df['bin1_start'] > close_diag_df['bin2_start']]

    assert (len(upper_triangle) + len(lower_triangle) == len(close_diag_df))

    upper_triangle_duplicates = upper_triangle[upper_triangle['count'] == upper_triangle['downleft']]
    upper_triangle_duplicates.loc[:,'downleft'] = 0

    print(upper_triangle_duplicates[['bin1_start', 'bin2_start'] + count_and_noise_columns])

    lower_triangle_duplicates = lower_triangle[lower_triangle['count'] == lower_triangle['upright']]
    lower_triangle_duplicates.loc[:,'upright'] = 0


    print(lower_triangle_duplicates[['bin1_start', 'bin2_start'] + count_and_noise_columns])

    noise_dataframe.loc[lower_triangle_duplicates.index, 'upright'] = lower_triangle_duplicates['upright']
    noise_dataframe.loc[upper_triangle_duplicates.index, 'downleft'] = upper_triangle_duplicates['downleft']

    #original = files.load_dataframe(noise_dataframe_path_1000)

    files.save_dataframe(noise_dataframe,noise_dataframe_path_1000)

if __name__ == "__main__":
    main()