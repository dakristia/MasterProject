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
    determine_breakpoints_for_counts_close_to_diag()

def determine_breakpoints_for_counts_close_to_diag():
    reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    reg_dataframe = files.load_dataframe(reg_dataframe_path_1000)

    reg_dataframe = reg_dataframe.sort_values(by=["modle_count"], ascending=False)

    result = reg_dataframe[(reg_dataframe['bin1_start'] - reg_dataframe['bin2_start']) > 1000].iloc[0]

    print(result)

    result = reg_dataframe[(reg_dataframe['bin1_start'] - reg_dataframe['bin2_start']) > 0].iloc[0]

    print(result)

if __name__ == "__main__":
    main()

