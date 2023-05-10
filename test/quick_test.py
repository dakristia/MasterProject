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

#@profile
def main():
    #correct_count_for_reg_pixels()
    #see_why_down_left_and_count_are_equal_2()
    #balanced_v_unbalanced()
    just_cooler_stuff()
    #cooler_fetch_example()


def cooler_fetch_example():
    import pandas as pd
    import numpy as np
    import cooler
    
    cooler_filepath = "../input/H1-hESC.mcool::/resolutions/1000"
    cooler_object = cooler.Cooler(cooler_filepath)

    dataframe_selector = cooler_object.matrix(as_pixels=True, join=True)
    matrix_selector = cooler_object.matrix(as_pixels=False, join=True)

    fetch_str = "chr1:1,000,000-10,000,000"

    dataframe = dataframe_selector.fetch(fetch_str)
    matrix = matrix_selector.fetch(fetch_str)


    print(dataframe)
    print(matrix)

def just_cooler_stuff():
    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"
    cooler_object = cooler.Cooler(cooler_file)

    dataframe_selector = cooler_object.matrix(as_pixels=True, join=False)
    matrix_selector = cooler_object.matrix(as_pixels = False,join=False)

    fetch_str = "chr1:1,000,000-1,500,000"

    print(type(dataframe_selector))

    dataframe = dataframe_selector.fetch(fetch_str)

    print(dataframe)

    matrix = matrix_selector.fetch(fetch_str)

    print(matrix)

    print(len(dataframe), matrix.shape[0] * matrix.shape[1])

    df_chrom1 = dataframe_selector.fetch("chr1")
    print(df_chrom1)
    print(df_chrom1.info(memory_usage="deep"))

    print(cooler_object.chromsizes)
    df_matrix1 = matrix_selector.fetch("chr1")
    #print(f"{round(sys.getsizeof(dataframe) / (1024 * 1024),2)}MB, {round(sys.getsizeof(matrix)/ (1024 * 1024),2)}MB")
    
    #print(dataframe.info(memory_usage="deep"))


def balanced_v_unbalanced():
    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"
    cre_file = "../input/H1-hESC.7group.PLSELSONLY.bed"

    cooler_object = cooler.Cooler(cooler_file)

    #balanced_selector = cooler_object.matrix(as_pixels=True, join=True, balance=True)

    #print(balanced_selector.fetch("chr1:1200000-1300000"))
    # unbalanced_selector = cooler_object.matrix(as_pixels=True, join=True, balanced=False)


    balanced_output_path = "./output/balVunbal/balanced_regelems.csv"
    unbalanced_output_path = "./output/balVunbal/unbalanced_regelems.csv"

    balanced_input_path = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"

    #balanced_df = files.load_dataframe(balanced_input_path)
    balanced_df = note_regulatory_interaction.note_prom_enh_inter_gwide(cooler_file,balanced=True,bed_file_path=cre_file,output_path=balanced_output_path)
    #print(balanced_df)
    #unbalanced_df = files.load_dataframe(unbalanced_output_path)
    unbalanced_df = note_regulatory_interaction.note_prom_enh_inter_gwide(cooler_file,balanced=False,bed_file_path=cre_file,output_path=unbalanced_output_path)

    print(balanced_df)
    print(unbalanced_df)

    # * Unbalanced dataframe has exact same counts as unbalanced counts in balanced_df

    #equal_df = balanced_df[balanced_df[["chrom","bin1_start","bin2_start"]] == unbalanced_df[["chrom","bin1_start","bin2_start"]]]

    # result_df = pd.merge(balanced_df, unbalanced_df, how='inner', on=["chrom","bin1_start","bin2_start"])
    # result_df = result_df.fillna(0)

    print(np.max(balanced_df['valid_count']))

    plt.scatter(balanced_df['modle_count'],balanced_df['valid_count'])
    plt.savefig('./output/balVunbal/correlation.png')
    plot_utils.show_plot('./output/balVunbal/correlation.png')
    files.viewplot()


def diagonal_contact_info():
    reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    reg_dataframe = files.load_dataframe(reg_dataframe_path_1000)

    # * Get center counts
    reg_dataframe = reg_dataframe[reg_dataframe["bin1_start"] - reg_dataframe["bin2_start"] == 0]
    print(reg_dataframe["modle_count"].describe())


def see_why_down_left_and_count_are_equal_2():
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
    print("percentile",np.percentile(distance_array,25))
    print("percentile",np.percentile(distance_array,50))
    print("percentile",np.percentile(distance_array,75))
    print("percentile",np.percentile(distance_array,100))

def correct_count_for_reg_pixels():
    reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    reg_dataframe = files.load_dataframe(reg_dataframe_path_1000)

    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"
    cooler_object = cooler.Cooler(cooler_file)
    selector = cooler_object.matrix(balance=False)

    correct_counts = 0
    incorrect_counts = 0
    surrounding_count_correct = 0

    resolution = 1000
    print(reg_dataframe.columns)
    for row in reg_dataframe.itertuples():

        df_index = int(row[0])
        bin1_start = int(row[8])
        bin1_end = int(row[9])
        bin2_start = int(row[10])
        bin2_end = int(row[11])
        count = int(row[12])
        chrom = row[1]
        #print(bin1_start)

        if bin1_start == bin2_start: continue
        if bin1_start > bin2_start:
            bin1_start, bin2_start = bin2_start, bin1_start
            bin1_end, bin2_end = bin2_end, bin1_end

        assert bin1_start <= bin2_start


        surrounding_matrix = selector.fetch(f"{chrom}:{bin1_start - resolution}-{bin1_end + resolution}",f"{chrom}:{bin2_start - resolution}-{bin2_end + resolution}")
        #print(df)
        matrix_count = surrounding_matrix[1][1]


        #print(count, matrix_count)
        if count == matrix_count:
            correct_counts += 1
        else:
            incorrect_counts += 1
            surround = 0
            for sur_row in surrounding_matrix:
                for sur_val in sur_row:
                    if sur_val == count:
                        surround = 1
            surrounding_count_correct += surround
    print("Reg pixel stats:")
    print("Correct counts:", correct_counts)
    print("Wrong counts:", incorrect_counts)
    print("Wrong, but correct in neighbor: ", surrounding_count_correct)



def see_why_down_left_and_count_are_equal():
    noise_dataframe_path_1000 = "./output/noise_dataframe_h1_1000.csv"
    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)
    reg_dataframe_path_1000 = "../input/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    reg_dataframe = files.load_dataframe(reg_dataframe_path_1000)

    cooler_file = "../input/H1-hESC.mcool::/resolutions/1000"

    over_800_noise_df = noise_dataframe[noise_dataframe['count'] > 800]

    pixels_over_800 = over_800_noise_df[["bin1_start","bin1_end","bin2_start","bin2_end","count", "chrom"]]
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
    for pixels in pixels_over_800.itertuples():
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


        matrix = selector.fetch(f"{chrom}:{bin1_start}-{bin1_end}",f"{chrom}:{bin2_start}-{bin2_end}")
        matrix_flip = selector.fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")
        df = df_selector.fetch(f"{chrom}:{bin1_start}-{bin1_end}",f"{chrom}:{bin2_start}-{bin2_end}")
        df_flip = df_selector.fetch(f"{chrom}:{bin2_start}-{bin2_end}",f"{chrom}:{bin1_start}-{bin1_end}")
        #print(df)
        matrix_count = matrix[0][0]
        #print(count, matrix_count)
        surrounding_matrix = selector.fetch(f"{chrom}:{bin1_start - resolution}-{bin1_end + resolution}",f"{chrom}:{bin2_start - resolution}-{bin2_end + resolution}")
        surrounding_matrix_flip = selector.fetch(f"{chrom}:{bin2_start - resolution}-{bin2_end + resolution}",f"{chrom}:{bin1_start - resolution}-{bin1_end + resolution}")
        # try:
        #     assert count == matrix_count
        # except AssertionError:
        #     print("ASSERTIONERROR",)
        #     print("count v matrix_count",count,matrix_count)
        #     print(pixels)
        #     print(df)
        #     print(df_flip)
        #     print(noise_dataframe.iloc[df_index])
        #     print(matrix)
        #     print(matrix_flip)
        #     print(surrounding_matrix)
        #     print(surrounding_matrix_flip)
        #     print()
        #     surrounding_equal = 0
        #     for sur_row in surrounding_matrix:
        #         for sur_val in sur_row:
        #             # print()
        #             # print(sur_val)
        #             # print(count)
        #             if sur_val == count:
        #                 surrounding_equal += 1

        #     assert surrounding_equal > 0
        #     # * If we get here, a pixel has been awarded the counts of its neighbor.
        #     # * This should be fixed elsewhere, and this test will proceed
        #     count = matrix_count
        # # * Check downleft

        # dl1_start = bin1_start + resolution
        # dl1_end = bin1_end + resolution
        # dl2_start = bin2_start - resolution
        # dl2_end = bin2_end - resolution


        # matrix = selector.fetch(f"{chrom}:{dl1_start}-{dl1_end}",f"{chrom}:{dl2_start}-{dl2_end}")
        # df = df_selector.fetch(f"{chrom}:{dl1_start}-{dl1_end}",f"{chrom}:{dl2_start}-{dl2_end}")
        # #print(df)
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
    #print("lowest count that has equal noise:", lowest_equal_count)
        #target_neighbor = 



def test_heatmap_and_overlay():
    import numpy as np
    import matplotlib.pyplot as plt # type:ignore

    # generate random data for the heatmap
    data = np.random.rand(10, 10)

    resolution = 5000

    cooler_object1 = cooler.Cooler(f"../input/H1-hESC.mcool::/resolutions/{resolution}")

    selector1 = cooler_object1.matrix(balance=False)
    selector2 = cooler_object1.matrix(balance=False,as_pixels=True,join=True)

    chrom = "chr19"

    bins = 9
    offset = 5000 * bins
    start = 1_000_000
    end = start + offset

    matrix = selector1.fetch((chrom,start,end))
    df = selector2.fetch((chrom,start,end))

    re_pair_data_path = f"../input/dataframes/H1_binsize{resolution}/genome_wide/{resolution}_H1-hESC.7group.csv"

    re_dataframe = files.load_dataframe(re_pair_data_path).query(f"chrom == chrom")

    query = re_dataframe.query(f'prom_start >= {start} and enh_start >= {start} and prom_end <= {end} and enh_end <= {end}')

    print(query)
    sub_matrix = matrix[:9,:9]
    #sub_matrix2 = matrix2[:9,:9]
    print(df)
    print(sub_matrix)
    #print(sub_matrix.dtype)

    np_matrix = np.asarray(sub_matrix)
    data = np_matrix

    # set the colormap to the one typically associated with Hi-C data
    cmap = 'YlOrRd'

    # create a figure and axis object
    fig, ax = plt.subplots()

    # create the heatmap with the 'none' interpolation option to hide the lines between cells
    heatmap = ax.imshow(data, cmap=cmap, interpolation='none')


    from matplotlib.colors import LinearSegmentedColormap # type:ignore
    # create a blue color map with positive values
    cmap2 = LinearSegmentedColormap.from_list("", ["lightblue", "blue"]) 
    #cmapWhite = LinearSegmentedColormap.from_list("", ["white", "white"]) 

    data2 = np.zeros((9, 9))
    data2[6, 6] = 6192
    data2[0, 6] = 213
    #print(data2)



    # set the x and y ticks to empty lists to hide the tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    # remove the borders and spines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # add a colorbar to the figure
    cbar = fig.colorbar(heatmap)

    from matplotlib.colors import LogNorm # type:ignore
    heatmap = ax.imshow(data, cmap=cmap, norm=LogNorm())

    # heatmapWhite = ax.imshow(data2, cmap=cmapWhite, norm=LogNorm())
    # heatmapWhite.set_alpha(np.where(data2 == 0, 0, 1))

    heatmap2 = ax.imshow(data2, cmap=cmap2, norm=LogNorm())
    heatmap2.set_alpha(np.where(data2 == 0, 0, 1.0))

    # show the plot
    plt.show()
    plt.savefig("Heatmap_TEST.png")
    plot_utils.show_plot("Heatmap_TEST.png")

def test_fetch_and_query():
    cooler_object1 = cooler.Cooler("../input/H1-hESC.mcool::/resolutions/5000")
    cooler_object2 = cooler.Cooler("../input/H1-hESC.mcool::/resolutions/1000")

    selector1 = cooler_object1.matrix(balance=False, as_pixels=True, join=True)
    selector2 = cooler_object2.matrix(balance=False, as_pixels=True, join=True)  

    offset = 10_000

    df1 = selector1.fetch(("chr1",1_000_000 + offset,1_010_000 + offset))
    df2 = selector2.fetch(("chr1",1_000_000 + offset,1_010_000 + offset))

    print(df1.query(f'start1 < 1_015_000 and start2 >= 1_015_000'))
    print(df2.query(f'start1 < 1_015_000 and start2 >= 1_015_000'))

    print(np.sum(df2.query(f'start1 < 1_015_000 and start2 >= 1_015_000')['count']))

def old_test():
    dataframe = pd.read_csv("../input/H1-hESC.7group.PLSELSONLY.bed",delim_whitespace=True,header=None,names=Constants.DATAFRAME_COLUMNS_BED)

    #dataframe = dataframe[dataframe["chrom"] == "chr1"]
    pls_dataframe = dataframe[dataframe["type"].str.contains('PLS')]
    els_dataframe = dataframe[dataframe["type"].str.contains('ELS')]

    chrom_names = np.unique(np.array(dataframe["chrom"]))

    count = 0

    for name in chrom_names:
        chrom_dataframe = dataframe[dataframe["chrom"] == name]
        pls_dataframe = dataframe[dataframe["type"].str.contains('PLS')]
        els_dataframe = dataframe[dataframe["type"].str.contains('ELS')]

        # start_column = np.array(chrom_dataframe["chromStart"])
        # end_column = np.array(chrom_dataframe["chromEnd"])
        print("number of promoters:", len(pls_dataframe))
        print("number of enhancers:", len(els_dataframe))

        print()
        for e in pls_dataframe.itertuples():
            for y in els_dataframe.itertuples():
                if abs(e[2] - y[3]) > 2_999_900: break
                if abs(e[3] - y[2]) > 2_999_900: break
                if abs(e[2] - y[3]) < 50 or abs(e[3] - y[2]) < 50 : count += 1
                #count += 1
                # if abs(e - y) <= 100:
                #     print(y, e)

        print(f"count after {name}: {count}")
        exit()
    print(count)
    exit()

    sizes = np.array([])

    for row in chrom_dataframe.itertuples():
        chromStart = row[2]
        chromEnd = row[3]
        size = abs( chromStart - chromEnd)
        #if chromStart % 100 == 50 and chromEnd % 100 == 50: print(row)
        #if size == 150: print(row)
        sizes = np.append(sizes,size )

    print(np.min(sizes))
    print(np.max(sizes))
    print(np.mean(sizes))
    print(np.median(sizes))

    print(np.count_nonzero())

    print(len(sizes))




if __name__ == "__main__":

    main()
    