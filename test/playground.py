import pandas as pd
import numpy as np
import cooler
from line_profiler import LineProfiler # type:ignore

import sys

sys.path.append('../src/regelem')
import plot_utils
import files
import Constants
import note_regulatory_interaction

def main():
    my_god()

def my_god():
    import matplotlib.pyplot as plt #type:ignore
    import numpy as np

    # Generate x-axis data
    x = np.linspace(0, 10, 1000)

    # Generate y-axis data for the descending line
    y1 = 100 - x * 10

    # Create the plot
    plt.plot(x, y1, label='Line 1')

    # Set the x-axis label and limits
    plt.xlabel('Time')
    plt.xlim([0, 10])

    # Set the y-axis label and limits
    plt.ylabel('Value')
    plt.ylim([0, 100])

    # Add a legend and show the plot
    plt.legend()
    plt.show()
    plt.savefig('stresslmao.png')
    plot_utils.show_plot('stresslmao.png')

def see_why_down_left_and_count_are_equal():
    noise_dataframe_path_1000 = "./output/noise_dataframe_h1_1000.csv"
    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)
    reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    reg_dataframe = files.load_dataframe(reg_dataframe_path_1000)

    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"

    over_587_noise_df = noise_dataframe[noise_dataframe['count'] > 587]

    pixels_over_587 = over_587_noise_df[["bin1_start","bin1_end","bin2_start","bin2_end","count", "chrom"]]
    reduced_noise_dataframe = noise_dataframe[["bin1_start","bin1_end","bin2_start","bin2_end","count", "chrom"]]
    reduced_noise_dataframe = reduced_noise_dataframe.sort_values(by=["count"])
    reduced_noise_dataframe = reduced_noise_dataframe.iloc[::-1]

    cooler_object = cooler.Cooler(cooler_file)
    selector = cooler_object.matrix(balance=False)
    df_selector = cooler_object.matrix(balance=False,as_pixels=True, join=True)

    n_equal = 0
    n_diff = 0
    lowest_equal_count = False

    resolution = 1000

    distance_array = np.array([])

    for pixels in pixels_over_587.itertuples():
        df_index = int(pixels[0])
        bin1_start = int(pixels[1])
        bin1_end = int(pixels[2])
        bin2_start = int(pixels[3])
        bin2_end = int(pixels[4])
        count = int(pixels[5])
        chrom = pixels[6]

        # * Skip pixels on diagonal
        if bin1_start == bin2_start: continue

        # * Look at pixels above diagonal
        if bin1_start > bin2_start:
            bin1_start, bin2_start = bin2_start, bin1_start
            bin1_end, bin2_end = bin2_end, bin1_end

        # * Ensure we did it right
        assert bin1_start <= bin2_start

        distance = bin2_start - bin1_start
        distance_array = np.append(distance_array, distance)

        # * Fetch 3x3 matrix surrounding 
        surrounding_matrix = selector.fetch(f"{chrom}:{bin1_start - resolution}-{bin1_end + resolution}",f"{chrom}:{bin2_start - resolution}-{bin2_end + resolution}")

        center_count = surrounding_matrix[1][1]
        matrix_downleft_noise = surrounding_matrix[2][0]
    
        #print(count, matrix_count)
        if center_count == matrix_downleft_noise:
            if not lowest_equal_count: lowest_equal_count = center_count
            elif lowest_equal_count > center_count: lowest_equal_count = center_count
            n_equal += 1
        else:
            n_diff += 1
            
    
    print(f"equal: {n_equal}")
    print(f"diff: {n_diff}")
    print(f"equal ratio: {n_equal / (n_equal + n_diff)}")
    print("total pixels checked:", n_equal + n_diff)

    #numpy.mean(), numpy.std(), numpy.min(), numpy.max(),  numpy.percentile()
    print()
    print("distance stats:",)
    print("mean",np.mean(distance_array))
    print("std",np.std(distance_array))
    print("min",np.min(distance_array))
    print("max",np.max(distance_array))
    print("percentile",np.percentile(distance_array))

if __name__ == "__main__":
    main()