import pandas as pd
import numpy as np
import cooler
from line_profiler import LineProfiler # type:ignore


import sys
sys.path.append('../src/regelem')
import plot_utils
import files
import Constants

#@profile
def main():

    test_heatmap_and_overlay()


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
    