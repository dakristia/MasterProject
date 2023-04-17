import pandas as pd
import numpy as np
import os
import sys
import time
import cooler

from line_profiler import LineProfiler #type:ignore

sys.path.append('../src/regelem')
import noise_analysis as na
import files


h1_dataframe = "../input/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
cooler_file = "../input/H1-hESC.mcool::/resolutions/5000"

def main():
    check_symmetry_of_cool(cooler_file)
    #validate_dataframe(h1_dataframe,cooler_file)

def check_symmetry_of_cool(path_to_cooler):
    cooler_object = cooler.Cooler(path_to_cooler)
    cooler_selector = cooler_object.matrix(balance=False,as_pixels=True,join=True)

    cooler_resolution = int(cooler_object.binsize)

    bin_counter = 0

    chrom = "chr19"

    for i in range(100):
        for y in range(100):
            bin1_start = i * cooler_resolution + 1_000_000
            bin1_end = bin1_start + cooler_resolution

            bin2_start = y * cooler_resolution + 1_000_000
            bin2_end = bin2_start + cooler_resolution

            dataframe = cooler_object.matrix(balance=False,as_pixels=True,join=True).fetch(f"{chrom}:{bin1_start}-{bin1_end}",f"{chrom}:{bin2_start}-{bin2_end}")
            dataframe2 = cooler_object.matrix(balance=False,as_pixels=True,join=True).fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")

            if len(dataframe):
                matrix = cooler_object.matrix(balance=False).fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")

                print()
                print()

                print()


                print(dataframe)

                print(dataframe2)

                print(matrix)
                print(dataframe.equals(dataframe2))


def validate_dataframe(path_to_dataframe : str, path_to_cooler : str):
    reg_dataframe = files.load_dataframe(path_to_dataframe)

    cooler_object = cooler.Cooler(path_to_cooler)
    cooler_selector = cooler_object.matrix(balance=False,as_pixels=True,join=True)

    cooler_resolution = int(cooler_object.binsize)


    number_of_bins_not_in_cooler = 0
    number_of_mismatches_bins = 0
    number_of_matched_bins = 0
    number_of_multiple_matches = 0


    for row in reg_dataframe.itertuples():
        
        chrom_index = reg_dataframe.columns.get_loc("chrom") + 1
        prom_name_index = reg_dataframe.columns.get_loc("prom_name") + 1
        enh_name_index = reg_dataframe.columns.get_loc("enh_name") + 1
        bin1_start_index = reg_dataframe.columns.get_loc("bin1_start") + 1
        bin1_end_index = reg_dataframe.columns.get_loc("bin1_end") + 1
        bin2_start_index = reg_dataframe.columns.get_loc("bin2_start") + 1
        bin2_end_index = reg_dataframe.columns.get_loc("bin2_end") + 1
        modle_count_index = reg_dataframe.columns.get_loc("modle_count") + 1

        chrom = row[chrom_index]
        prom_name = row[prom_name_index]
        enh_name = row[enh_name_index]
        bin1_start = row[bin1_start_index]
        bin1_end = row[bin1_end_index]
        bin2_start = row[bin2_start_index]
        bin2_end = row[bin2_end_index]
        count = row[modle_count_index]

        resolution = int(bin1_end - bin1_start)


        if resolution != cooler_resolution:
            print(f"ERROR: Resolution of dataframe {resolution} does not match resolution of cooler file {cooler_resolution}")
            return

        if resolution != int(bin2_end - bin2_start):
            print(f"ERROR: Resolution of bin 1 {resolution} does not match resolution of bin 2 {bin2_end - bin2_start}")
            return

        x = bin1_start // resolution
        y = bin2_start // resolution

        #dataframe = cooler_object.matrix(balance=False,as_pixels=True,join=True).fetch(f"{chrom}:{bin1_start}-{bin1_end}",f"{chrom}:{bin2_start}-{bin2_end}")
        dataframe = cooler_object.matrix(balance=False,as_pixels=True,join=True).fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")

        matrix = cooler_object.matrix(balance=False,as_pixels=True,join=True).fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")

        # print(dataframe)
        # print(row)
        # print(dataframe2)

        # if not dataframe.equals(dataframe2):
        #     print(row)
        #     print(dataframe)
        #     print(dataframe2)

        # continue
        # print(dataframe)
        # print(pd.DataFrame(row))
        # if len(dataframe): exit()

        last_mismatch = number_of_mismatches_bins

        if len(dataframe) == 1:
            cooler_count = np.array(dataframe["count"])[0]
            if cooler_count == count:
                number_of_matched_bins += 1
            else:
                number_of_mismatches_bins += 1
        elif len(dataframe) == 0:
            number_of_bins_not_in_cooler += 1
        else:
            cooler_counts = np.array(dataframe["count"])
            number_of_multiple_matches += 1
            
            all_match = True
            for n in cooler_counts:
                if n != count:
                    number_of_mismatches_bins += 1
                    all_match = False
                    break

            if all_match: number_of_matched_bins += 1

        if last_mismatch != number_of_mismatches_bins:
            print("New mistmatch!")


    dataframe_size = len(reg_dataframe)

    print(f"Length of dataframe: {dataframe_size}")
    print()
    print(f"Matched: {number_of_matched_bins} {number_of_matched_bins/dataframe_size}")
    print(f"Mismatched: {number_of_matched_bins} {number_of_mismatches_bins/dataframe_size}")
    print(f"Missing: {number_of_bins_not_in_cooler} {number_of_bins_not_in_cooler/dataframe_size}")
    print(f"Multiple match: {number_of_multiple_matches} {number_of_multiple_matches/dataframe_size}")


if __name__ == "__main__":
    main()