import pandas as pd
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt #type:ignore

#from test_noise_analysis import plot_noise
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

    output_csv_path = "./output/noise_dataframe_h1_1000_removeDuplicateCloseToDiagonal.csv"
    #noise_dataframe = na.register_noise(reg_dataframe_path_1000,cooler_path,res,output_csv_path)

    noise_dataframe = files.load_dataframe(noise_dataframe_path_1000)
    lesser_res_dataframe = files.load_dataframe(reg_dataframe_path)

    plot_noise(noise_dataframe,lesser_res_dataframe,output_plot_path_1000)
    #more_noise_plotting(noise_dataframe)
    # #profile_funct()
    # exit()

def plot_noise(noise_dataframe : pd.DataFrame, lesser_res_dataframe : pd.DataFrame, output_path = False):
    ""


    noise_dataframe = pd.DataFrame(noise_dataframe).sort_values(by=["chrom","bin1_start","bin2_start"])

    # Remove duplicate pixels. This will remove some pairs, but we only care about specific pixels right now. 
    duplicates = noise_dataframe.loc[:,["chrom","bin1_start","bin1_end","bin2_start","bin2_end"]].duplicated(keep="first",)
    
    noise_dataframe = noise_dataframe[~duplicates].sort_values(by=["chrom","bin1_start","bin2_start"])

    chrom_array = np.array(noise_dataframe["chrom"])
    count_array = np.array(noise_dataframe["count"])
    noise_array = np.array(noise_dataframe["noise"])
    noise_count_ratio_array = noise_array / count_array

    adjacent_columns = ["down","downleft","left","upleft","up","upright","right","downright"]


    
    max_noise_array = np.array(noise_dataframe[adjacent_columns].max(axis=1))
    min_noise_array = np.array(noise_dataframe[adjacent_columns].min(axis=1))
    average_noise_array = np.array(noise_dataframe[adjacent_columns].mean(axis=1))
    std_noise_array = np.array(noise_dataframe[adjacent_columns].std(axis=1))
    adjacent_pixels_array = np.array(noise_dataframe['adjacent'])
    distance_array = np.array(abs(noise_dataframe['bin1_start'] - noise_dataframe['bin2_start']))


    ## * Sorted version of the arrays
    sorted_count_array, sorted_noise_array, sorted_std_array, sorted_distance_array = zip(*sorted(zip(count_array, noise_array, std_noise_array, distance_array)))
    sorted_count_array = np.array(sorted_count_array)
    sorted_noise_array = np.array(sorted_noise_array)
    sorted_std_array = np.array(sorted_std_array)

    average_sorted_noise_array = sorted_noise_array/adjacent_pixels_array

    # create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12,10))
    
    # create a set of unique labels
    unique_labels = set(chrom_array)
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap #type:ignore
    viridis = plt.cm.get_cmap('gist_rainbow', len(unique_labels))
    color_choices = viridis(np.linspace(0, 1, 24))
    color_dict = {label: color_choices[i] for i, label in enumerate(unique_labels)}
    colors = [color_dict[label] for label in chrom_array]
    
    # * Plot 1

    for label in unique_labels:
        indices = [i for i, x in enumerate(chrom_array) if x == label]
        uq_count_array = [count_array[i] for i in indices]
        uq_noise_array = [noise_array[i]/8 for i in indices]
        axs[0, 0].scatter(uq_count_array, uq_noise_array, alpha = 0.4, label=label ,color=color_dict[label], s = 20)

    axs[0, 0].plot(count_array, count_array, color="black")

    axs[0, 0].set_title('Count vs noise per pair')
    axs[0, 0].set_xlabel('Center frequency "count"')
    axs[0, 0].set_ylabel('Average surrounding\nfrequency "noise"')

    plt.rcParams['legend.fontsize'] = 'small'
    axs[0, 0].legend(loc='upper center', columnspacing=0.2, bbox_to_anchor=(0.22, 1), ncol=3, prop={'size': 6})

    ######### * Plot 2
    #axs[0, 1].scatter(sorted_count_array, average_sorted_noise_array, color="grey")

    # * FOURTH trend line
    # * Define the cubic function to fit
    def cubic(sorted_count_array, a, b, c, d):
        return a*sorted_count_array**3 + b*sorted_count_array**2 + c*sorted_count_array + d


    # * Define the quadratic function to fit. Unused.
    def quadratic(sorted_count_array, a, b, c):
        return a*sorted_count_array**2 + b*sorted_count_array + c

    from scipy.optimize import curve_fit
    # * Fit the quadratic function to the data
    # * Simply change this and y_fit further down to use quadratic instead if we want a different curve.

    import seaborn as sns #type:ignore
    from scipy.stats import gaussian_kde

    density_plot = False

    if density_plot:

        # * Claculate the density of the points
        vstacked_flatten_counts = np.vstack([sorted_count_array,average_sorted_noise_array])
        gaussian = gaussian_kde(vstacked_flatten_counts)(vstacked_flatten_counts)
        # * Sort the points by density, so that the densest points are plotted last
        idx = gaussian.argsort()
        x, y, gaussian = sorted_count_array[idx], average_sorted_noise_array[idx], gaussian[idx]

        count_noise_dataframe = pd.DataFrame({"df1":sorted_count_array,"df2":average_sorted_noise_array})

        sns.scatterplot(
            data=count_noise_dataframe,
            x = "df1",
            y = "df2",
            c=gaussian,
            cmap="viridis",
            ax=axs[0,1],
        )

    else:
        #* Create a scatterplot of the data
        axs[0,1].scatter(sorted_count_array, average_sorted_noise_array, color="grey")


    # * Set y axis to match previous plot
    axs[0, 1].set_ylim(axs[0, 0].get_ylim())



    # * Create a line plot of the fitted cubic function
    # * We plot a third-order polynomial to this data, as it seems to fit the best.
    popt, pcov = curve_fit(cubic, sorted_count_array, average_sorted_noise_array) 
    x_fit = np.linspace(sorted_count_array.min(), sorted_count_array.max(), 100)
    y_fit = cubic(x_fit, *popt)
    axs[0,1].plot(x_fit, y_fit, 'r-', label='Running trend')

    # * Calculate the standard deviation of the data
    #! Unused
    residuals = average_sorted_noise_array - cubic(sorted_count_array, *popt)
    stdev = np.std(residuals)

    # * This better explains the variation of the data than stdev. 
    std_between = np.std(np.concatenate((sorted_count_array, average_sorted_noise_array)))

    # * Add upper and lower standard deviation curves
    axs[0,1].plot(x_fit, y_fit + std_between, 'b--', label='Upper Std Dev')
    axs[0,1].plot(x_fit, y_fit - std_between, 'g--', label='Lower Std Dev')
    axs[0,1].legend()

    axs[0, 1].set_title('Count vs noise per pair')
    axs[0, 1].set_xlabel('Center frequency "count"')
    axs[0, 1].set_ylabel('Average surrounding\nfrequency "noise"')
    

    # * Plot 3

    #! Get standardeviation of noise for each pixel

    for label in unique_labels:
        indices = [i for i, x in enumerate(chrom_array) if x == label]
        uq_count_array = [count_array[i] for i in indices]
        uq_std_array = [std_noise_array[i] for i in indices]
        scatter = axs[1, 0].scatter(uq_count_array, uq_std_array, alpha = 0.4, label=label ,color=color_dict[label], s = 20)

    corr_coef = np.corrcoef(sorted_distance_array, average_sorted_noise_array)[0, 1]
    print("Correlation Coefficient:", corr_coef)


    axs[1,0].set_title('Count vs std noise per pair')
    axs[1,0].set_xlabel('Center frequency "count"')
    axs[1,0].set_ylabel('Standard deviation of\nsurrounding frequency "noise"')

    popt, pcov = curve_fit(cubic, sorted_count_array, sorted_std_array) 
    x_fit = np.linspace(sorted_count_array.min(), sorted_count_array.max(), 100)
    y_fit = cubic(x_fit, *popt)
    axs[1,0].plot(x_fit, y_fit, 'r-')

    # coeffs = np.polyfit(sorted_count_array, sorted_std_array, 2)
    # p = np.poly1d(coeffs)
    # axs[1,0].plot(sorted_count_array, p(sorted_count_array), color='red')

    #slope, intercept = np.polyfit(sorted_count_array, average_sorted_noise_array, 1)
    #axs[1,0].plot(sorted_count_array, slope * sorted_count_array + intercept, color='red')
    
    # * Reduce font size in legend
    plt.rcParams['legend.fontsize'] = 'small'

    # * Plot 4
    #sub_gs = axs[1, 1].subgridspec(2, 2)
    #axs[1, 1] = plt.subplots(1,2)
    #sub_axs = axs[1, 1].subplots(2)
    #from mpl_toolkits.axes_grid1.inset_locator import inset_axes
   # positions = [0.4, 0.8, 1.2, 1.6]
    #ax_inset = inset_axes(axs[1, 1], width="50%", height="50%", loc="lower left", borderpad=1)
    #sub_axs[0,0].boxplot([count_array, min_noise_array, max_noise_array, average_noise_array], labels=["Count","Min\nnoise","Max\nnoise","Average\nnoise"], positions=positions, widths=0.3)
    axs[1,1].boxplot([count_array, min_noise_array, max_noise_array, average_noise_array], labels=["Count","Min\nnoise","Max\nnoise","Average\nnoise"],  widths=0.3)
    
    #ax_inset.boxplot([count_array, min_noise_array, max_noise_array, average_noise_array], labels=["A","B","C","D"], positions=[p+1.6 for p in positions], widths=0.3)
    # #axs
    # axs[1,1].set_xticks([1, 2])
    # axs[1,1].set_xlabel('Scale')
    # axs[1,1].set_ylabel('Value')
    axs[1,1].set_title('Boxplots')

    # Set the scale of the second boxplot to logarithmic
    axs[1,1].set_yscale('symlog')

    # Adjust the layout to make room for the larger x-axis label
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    #axs[1, 1].scatter(sorted_distance_array, average_sorted_noise_array, color='orange')
    
    corr_coef = np.corrcoef(sorted_distance_array, average_sorted_noise_array)[0, 1]

    #print("Correlation Coefficient:", corr_coef)
    #axs[1, 1].scatter(sorted_std_array,average_sorted_noise_array, color='orange')
    #axs[1, 1].set_xlabel('Distance between bins')
    #axs[1, 1].set_ylabel('Average surrounding\nfrequency "noise"')
    #box = axs[1, 1].boxplot([count_array, min_noise_array, max_noise_array, average_noise_array], labels=["Count","Min\nnoise","Max\nnoise","Average\nnoise"])
    #axs[1, 1].set_xscale('log')
    #axs[1, 1].set_yscale('log')

    

    # * Tried to find min noise with this, but the box simply isnt there
    # ymin, ymax =  axs[1, 1].get_ylim()
    # axs[1, 1].set_ylim([min(ymin, 1e-6), ymax])

    # # add a horizontal line at y=0
    # axs[1, 1].axhline(y=1e-6, color='k', linestyle='--')

    # * Show the plot
    plt.show()
    if output_path: 
        plt.savefig(output_path)
        plot_utils.show_plot(output_path)

def more_noise_plotting(noise_dataframe : pd.DataFrame):

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

    # * Create a column 'is_duplicate' to mark duplicates
    noise_dataframe = noise_dataframe.assign(is_duplicate=noise_dataframe.loc[:, ["chrom", "bin1_start", "bin1_end", "bin2_start", "bin2_end"]].duplicated(keep="first"))

    # * Group by the columns of interest and count the duplicates
    duplicate_counts = noise_dataframe.groupby(["chrom","bin1_start", "bin1_end", "bin2_start", "bin2_end"]).agg({'is_duplicate': 'sum'}).reset_index()

    # * Rename the 'is_duplicate' column to 'duplicate_count'
    duplicate_counts = duplicate_counts.rename(columns={'is_duplicate': 'duplicate_count'})

    # * Remove duplicates from the original DataFrame
    noise_dataframe = noise_dataframe[~noise_dataframe['is_duplicate']].sort_values(by=["chrom", "bin1_start", "bin2_start"])

    # * Drop the 'is_duplicate' column from the original DataFrame
    noise_dataframe = noise_dataframe.drop(columns=['is_duplicate'])

    # * Merge the original DataFrame with the duplicate_counts DataFrame
    noise_dataframe = noise_dataframe.merge(duplicate_counts, on=["chrom", "bin1_start", "bin1_end", "bin2_start", "bin2_end"], how="left")

    # * Fill any NaN values in the 'duplicate_count' column with 0
    noise_dataframe['duplicate_count'] = noise_dataframe['duplicate_count'].fillna(0)

    # * Convert the 'duplicate_count' column to an integer type
    noise_dataframe['duplicate_count'] = noise_dataframe['duplicate_count'].astype(int)

    noise_dataframe['duplicate_count'] += 1
    noise_dataframe = noise_dataframe.rename(columns={'duplicate_count': 'npairs'})
    
    assert original_len == sum(noise_dataframe["npairs"])
    
    # * Get an overview of directional columns
    adjacent_columns = ["down","downleft","left","upleft","up","upright","right","downright"]
    
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

    # # * If bin1_start != bin2_start, and if [bin1_start,bin2_start], then [bin2_start,bin1_start] does not exist
    # original = non_diagonal_dataframe[["bin1_start","bin1_end","bin2_start","bin2_end"]]
    # swapped = original.copy()
    # swapped["bin1_start"], swapped["bin2_start"] = swapped["bin2_start"], swapped["bin1_start"]    
    # swapped["bin1_end"], swapped["bin2_end"] = swapped["bin2_end"], swapped["bin1_end"]    
    # # * Checks if mirrored exists
    # matching_rows = original.merge(swapped, how='inner')
    # assert len(matching_rows) == 0
    # # * Checks for duplicates
    # assert not (swapped["bin1_start"] == original["bin1_start"]).any()
    # assert not (swapped["bin1_end"] == original["bin1_end"]).any()

    breakoffs_array = [200,400,600,800,1000]
    prev_breakoff = -1

    breakoff_dataframes = np.empty(shape=len(breakoffs_array) + 1, dtype=object)

    for i,breakoff in enumerate(breakoffs_array):
        part_dataframe = non_diagonal_dataframe[non_diagonal_dataframe["count"] > prev_breakoff][non_diagonal_dataframe["count"] <= breakoff]
        breakoff_dataframes[i] = part_dataframe
        prev_breakoff = breakoff
    
    i += 1
    part_dataframe = non_diagonal_dataframe[non_diagonal_dataframe["count"] > prev_breakoff]
    breakoff_dataframes[i] = part_dataframe

    def plot_average_of_surround(dataframe_array : np.array,out_path : str, breakoffs : np.array):
        dataframe_len = len(dataframe_array)
        n_of_cols = 2
        n_of_rows = dataframe_len / n_of_cols
        if n_of_rows % 1 != 0: 
            n_of_rows = round(n_of_rows) + 1
            uneven = True
        else: 
            n_of_rows = int(n_of_rows)
            uneven = False
        
        print(n_of_rows)

        fig, axs = plt.subplots(nrows=n_of_rows, ncols=n_of_cols, figsize=(14, 10), dpi=300) 

        fig.subplots_adjust(wspace=0.3)

        cur_row = 0
        cur_col = 0

        bottom_range = 0
        upper_range = breakoffs[0]

    

        for index, df in enumerate(dataframe_array):
            if cur_col == n_of_cols: 
                cur_row += 1
                cur_col = 0

            mean_left = round(np.mean(df['left']),2)
            mean_upleft = round(np.mean(df['upleft']),2)
            mean_up = round(np.mean(df['up']),2)
            mean_upright = round(np.mean(df['upright']),2)
            mean_right = round(np.mean(df['right']),2)
            mean_downright = round(np.mean(df['downright']),2)
            mean_down = round(np.mean(df['down']),2)
            mean_downleft = round(np.mean(df['downleft']),2)
            mean_center = round(np.mean(df['count']),2)


            number_of_pixels = len(df)
            if upper_range: title_string = f"Counts: {bottom_range}-{upper_range}"
            else: title_string = f"Counts: {bottom_range} and above"
            title_string += f"\nNumber of pixels: {number_of_pixels}"


            heatmap_array = np.array([  [mean_upleft,mean_up,mean_upright],
                                        [mean_left,mean_center,mean_right],
                                        [mean_downleft,mean_down,mean_downright]])

            from matplotlib.colors import LogNorm #type:ignore
            im = axs[cur_row][cur_col].imshow(heatmap_array,cmap='magma', )#norm=LogNorm())
            axs[cur_row][cur_col].set_title(f'{title_string}')
            axs[cur_row][cur_col].set_xticks([])
            axs[cur_row][cur_col].set_yticks([])

            import matplotlib.patheffects as path_effects #type:ignore

            for i in range(heatmap_array.shape[0]):
                for j in range(heatmap_array.shape[1]):
                    text = axs[cur_row][cur_col].text(j, i, heatmap_array[i, j],
                        ha="center", va="center", color="w")
                    text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                                path_effects.Normal()])



            if index + 1 < len(breakoffs):
                bottom_range = breakoffs[index]
                upper_range = breakoffs[index+1]
            else:
                bottom_range = breakoffs[-1]
                upper_range = False

            cur_col +=1

        #axs[2, 1].set_visible(False)

        fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.6)
        fig.suptitle('Mean counts of center pixels and surrounding noise')
        plt.show()
        plt.savefig(out_path)

        plot_utils.show_plot(out_path)

    plot_average_of_surround(breakoff_dataframes, './output/plots/non_diagonal_data.png', breakoffs_array)

    # * Do the same for DIAGONAL data:

    prev_breakoff = -1

    breakoff_dataframes = np.empty(shape=len(breakoffs_array) + 1, dtype=object)

    for i,breakoff in enumerate(breakoffs_array):
        part_dataframe = diagonal_dataframe[diagonal_dataframe["count"] > prev_breakoff][diagonal_dataframe["count"] <= breakoff]
        breakoff_dataframes[i] = part_dataframe
        prev_breakoff = breakoff
    
    i += 1
    part_dataframe = diagonal_dataframe[diagonal_dataframe["count"] > prev_breakoff]
    breakoff_dataframes[i] = part_dataframe

    plot_average_of_surround(breakoff_dataframes,"./output/plots/diagonal_noise.png", breakoffs_array)

    #TODO: Correlation between number of pairs and noise?

    recollected_dataframe = pd.concat([diagonal_dataframe,non_diagonal_dataframe],ignore_index=True).sort_values(by=["chrom","bin1_start","bin2_start"])

    noise_array = np.array(recollected_dataframe["noise"])
    pair_array = np.array(recollected_dataframe['npairs'])

    plt.figure()

    noise_array = np.array(recollected_dataframe["noise"])
    pair_array = np.array(recollected_dataframe['npairs'])

    unique_pairs = np.unique(pair_array)  # get unique values of pair_array
    average_per_pair = np.empty(shape=unique_pairs.shape)
    for i,p in enumerate(unique_pairs):
        avg_noise = np.mean(noise_array[pair_array == p])
        average_per_pair[i] = avg_noise

    # * Plot number of pairs v average pers

    plt.plot(unique_pairs,average_per_pair,color='orange')
    plt.title('Average noise based on number of \npromoter/enhancer pers per pixel',fontsize=14)
    plt.xlabel('pairs per pixel',fontsize=12)
    plt.ylabel('noise',fontsize=12)

    plt.savefig('./output/plots/noise_v_pairs.png')
    plot_utils.show_plot('./output/plots/noise_v_pairs.png')



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