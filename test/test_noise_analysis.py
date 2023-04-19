import pandas as pd
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt #type:ignore

from line_profiler import LineProfiler #type:ignore

sys.path.append('../src/regelem')
import noise_analysis as na
import files
import plot_utils

cooler_path = "../input/H1-hESC.mcool::/resolutions/1000"
output_csv_path = "./output/noise_dataframe_1000.csv"
output_plot_path = "./output/plots/1000_noise_plot.png"
output_plot_path_1000 = "./output/plots/1000_noise_plot.png"

noise_dataframe_path_5000 = "./output/noise_dataframe_h1_5000.csv"
noise_dataframe_path_1000 = "./output/noise_dataframe_h1_1000.csv"
reg_dataframe_path = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"

res = 1000
# h1_dataframe = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
# cooler_file = "../input/H1-hESC.mcool::/resolutions/5000"


def main():
    #na.generate_bias("../input/H1-hESC.mcool::/resolutions/5000",output_path="bias_maybe.tsv")

    #noise_dataframe = na.register_noise(reg_dataframe_path_1000,cooler_path,res,output_csv_path)
    #exit()
    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)
    #lesse_res_dataframe = files.load_dataframe(reg_dataframe_path)

    look_more_at_low_frequencies(noise_dataframe)
    #na.plot_noise(noise_dataframe,lesse_res_dataframe,output_plot_path_1000)
    # #profile_funct()
    # exit()

def look_more_at_low_frequencies(noise_dataframe : pd.DataFrame):

    # * Divide dataframe into 10 percent segments. 
    # * First has the rows where 'count' is between 0 and 10 percent of the max value in 'count' column.
    # * Second has the rows between 10 and 20, third between 20 and 30 etc.

    def unmirror_diagonal(row):
        """Move pixel from the lower triangle of a matrix' diagonal to the upper triangle

        Args:
            row (_type_): DataFrame row

        Returns:
            _type_: _description_
        """
        if row['bin1_start'] > row['bin2_start'] and row['bin1_end'] > row['bin2_end']:
            row['bin1_start'], row['bin2_start'] = row['bin2_start'], row['bin1_start']
            row['bin1_end'], row['bin2_end'] = row['bin2_end'], row['bin1_end']
            row['left'],row['up'] = row['up'],row['left']
            row['downleft'],row['upright'] = row['upright'],row['downleft']
            row['down'],row['right'] = row['right'],row['down']
        return row


    def split_diagonal_from_rest(dataframe : pd.DataFrame) -> tuple:
        """Split dataframe into two new dataframes. One with every pixel on the diagonal, and one with every pixel not on diagonal

        Args:
            dataframe (pd.DataFrame): original dataframe

        Returns:
            tuple: Tuple containing diagonal dataframe on index 0, and non-diagonal dataframe on index 1
        """
        diagonal_dataframe = dataframe[(dataframe["bin1_start"] == dataframe["bin2_start"]) & (dataframe["bin1_end"] == dataframe["bin2_end"])]

        boolean_index = ~dataframe.isin(diagonal_dataframe).all(axis=1)

        non_diagonal_dataframe = dataframe[boolean_index]
        
        return diagonal_dataframe, non_diagonal_dataframe

    original_len = len(noise_dataframe)

    # * Remove duplicates based on bins only. We dont care about individual pairs (for now)'
    # TODO: Add column with number of promoter-enhancer pairs per pixel

    # * These columns are not needed for now
    noise_dataframe = noise_dataframe.drop(["prom_name","enh_name"],axis=1)

    # Create a column 'is_duplicate' to mark duplicates
    noise_dataframe = noise_dataframe.assign(is_duplicate=noise_dataframe.loc[:, ["chrom", "bin1_start", "bin1_end", "bin2_start", "bin2_end"]].duplicated(keep="first"))

    # Group by the columns of interest and count the duplicates
    duplicate_counts = noise_dataframe.groupby(["chrom","bin1_start", "bin1_end", "bin2_start", "bin2_end"]).agg({'is_duplicate': 'sum'}).reset_index()

    # Rename the 'is_duplicate' column to 'duplicate_count'
    duplicate_counts = duplicate_counts.rename(columns={'is_duplicate': 'duplicate_count'})

    # Remove duplicates from the original DataFrame
    noise_dataframe = noise_dataframe[~noise_dataframe['is_duplicate']].sort_values(by=["chrom", "bin1_start", "bin2_start"])

    # Drop the 'is_duplicate' column from the original DataFrame
    noise_dataframe = noise_dataframe.drop(columns=['is_duplicate'])

    # Merge the original DataFrame with the duplicate_counts DataFrame
    noise_dataframe = noise_dataframe.merge(duplicate_counts, on=["chrom", "bin1_start", "bin1_end", "bin2_start", "bin2_end"], how="left")

    # Fill any NaN values in the 'duplicate_count' column with 0
    noise_dataframe['duplicate_count'] = noise_dataframe['duplicate_count'].fillna(0)

    # Convert the 'duplicate_count' column to an integer type
    noise_dataframe['duplicate_count'] = noise_dataframe['duplicate_count'].astype(int)

    noise_dataframe['duplicate_count'] += 1
    noise_dataframe = noise_dataframe.rename(columns={'duplicate_count': 'npairs'})
    

    assert original_len == sum(noise_dataframe["npairs"])
    
    # * Get an overview of directional columns
    adjacent_columns = ["left","upleft","up","upright","right","downright","down","downleft"]
    
    # * lower_columns pixels are swapped with higher_columns when mirrored on the diagonal
    # * Variables unused, but serve as a note of sorts
    lower_columns = ["left","downleft","down"]
    higher_columns = ["up","upright","right"]

    # * Make sure all pixels are on the same side of the diagonal. 
    # * After this, numbers in higher_columns move away from the diagonal,
    # * while lower_columns move towards the diagonal. upleft, and downright dont change. 
    noise_dataframe = noise_dataframe.apply(unmirror_diagonal, axis=1)

    # * Split into diagonal and non-diagonal dataframes
    diagonal_dataframe, non_diagonal_dataframe = split_diagonal_from_rest(noise_dataframe)

    # * If bin1_start != bin2_start, and if [bin1_start,bin2_start], then [bin2_start,bin1_start] does not exist
    original = non_diagonal_dataframe[["bin1_start","bin1_end","bin2_start","bin2_end"]]
    swapped = original.copy()
    swapped["bin1_start"], swapped["bin2_start"] = swapped["bin2_start"], swapped["bin1_start"]    
    swapped["bin1_end"], swapped["bin2_end"] = swapped["bin2_end"], swapped["bin1_end"]    
    # * Checks if mirrored exists
    matching_rows = original.merge(swapped, how='inner')
    assert len(matching_rows) == 0
    # * Checks for duplicates
    assert not (swapped["bin1_start"] == original["bin1_start"]).any()
    assert not (swapped["bin1_end"] == original["bin1_end"]).any()

    # * Non-diagonal dataframe first
    array_ten_percent_dataframes = np.empty(shape=(10), dtype=pd.DataFrame)

    count_array = non_diagonal_dataframe["count"]

    highest_count = np.max(count_array)

    for i in range(0, 10):
        percentage_floor = i * 0.10 * highest_count
        percentage_roof = percentage_floor + 0.10 * highest_count
        ten_percent_dataframe = non_diagonal_dataframe[non_diagonal_dataframe["count"] > percentage_floor]
        ten_percent_dataframe = ten_percent_dataframe[ten_percent_dataframe["count"] <= percentage_roof]
        array_ten_percent_dataframes[i] = ten_percent_dataframe

    
    def plot_ten_percent_dataframes(array_ten_percent_dataframes, subtitle : str, output_path : str):
        # * Determine how many rows and columns that are needed to display the data
        dataframe_len = len(array_ten_percent_dataframes)
        n_of_cols = 2
        n_of_rows = dataframe_len / n_of_cols
        if n_of_rows % 1 != 0: n_of_rows = round(n_of_rows) + 1
        else: n_of_rows = int(n_of_rows)

        fig, axs = plt.subplots(nrows=n_of_rows, ncols=n_of_cols, figsize=(20, 20))

        cur_row = 0
        cur_col = 0

        describe_columns = ["count","noise","npairs"] + adjacent_columns
        
        low_percent = 0
        high_percent = 10
        low_count = 0
        high_count = 0

         # * loop through each dataframe and create a boxplot of its data
        for i, df in enumerate(array_ten_percent_dataframes):
            reduced_df = df[describe_columns]
            if cur_col == n_of_cols: 
                cur_row += 1
                cur_col = 0

            low_count = np.min(df["count"])
            high_count = np.max(df["count"])
            title_string = f"{low_percent}-{high_percent}% Counts: {low_count}-{high_count}"
            
            ax = reduced_df.plot(kind='box', ax=axs[cur_row][cur_col], title=f'{title_string}', logy=True)
            cur_col += 1
            low_percent += 10
            high_percent += 10

        fig.suptitle(subtitle)
        plt.savefig(output_path)
        plot_utils.show_plot(output_path)

    plot_ten_percent_dataframes(array_ten_percent_dataframes,'Non-diagonal data',"./output/non_diagonal_data.png")


    # * Do the same for DIAGONAL data:
    array_ten_percent_dataframes = np.empty(shape=(10), dtype=pd.DataFrame)

    count_array = diagonal_dataframe["count"]

    highest_count = np.max(count_array)

    for i in range(0, 10):
        percentage_floor = i * 0.10 * highest_count
        percentage_roof = percentage_floor + 0.10 * highest_count
        ten_percent_dataframe = diagonal_dataframe[diagonal_dataframe["count"] > percentage_floor]
        ten_percent_dataframe = ten_percent_dataframe[ten_percent_dataframe["count"] <= percentage_roof]
        array_ten_percent_dataframes[i] = ten_percent_dataframe

    plot_ten_percent_dataframes(array_ten_percent_dataframes,"Diagonal data", "./output/diagonal_noise.png")


    #TODO: Correlation between number of pairs and noise?


def profile_funct():
    # Create a new LineProfiler object
    profiler = LineProfiler()

    # Add the function to be profiled to the profiler
    profiler.add_function(na.register_noise)

    # Run the profiler on the function with arguments
    profiler.run('na.funct(reg_dataframe_path,cooler_path,res)')

    # Print the results
    profiler.print_stats()

if __name__ == "__main__":
    main()