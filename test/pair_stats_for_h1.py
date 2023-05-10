import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cooler
from line_profiler import LineProfiler # type:ignore


import sys
sys.path.append('../src/regelem')
import plot_utils
import files
import Constants
import note_regulatory_interaction


def main():
    input_data = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    df = files.load_dataframe(input_data)

    print(df)

    mean_count_per_bin = np.mean(df['modle_count'])
    print(mean_count_per_bin)

    


# However, the highest resolution available to us in the H1 data set we retrieved is 1000bp. 
# At this resolution, we see an average of (TODO: Insert average) counts per regulatory bin, 
# and a maximum of (TODO: Insert maximum number of pairs per bin) enhancer-promoter pairs per bin. 
# The frequency of multi-pair bins is (TODO: Insert multi-pair bins). 
# This does mean that (TODO: percentage of single-pair bins) percent of bins contain only one regulatory pair. 
# However, at this resolution, it is possible for a bin to contain several other elements that are not 
# enhancers and promoters as well. When we are looking at Micro-C data, it may be unclear if a single-pair bin 
# shows the interaction between the pair in question, or some other elements present in the bin.


def pixels_on_diagonal(chrom_size : int, diagonal_width: int):
    n_of_pixels = 0

    for i in range(chrom_size-diagonal_width):
            n_of_pixels += 2 * diagonal_width + 1
    for i in range(diagonal_width):
            n_of_pixels += 1 + i * 2
    
    return n_of_pixels


if __name__ == "__main__":
    main()
