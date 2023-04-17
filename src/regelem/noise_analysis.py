import pandas as pd
import numpy as np
import cooler
import Constants
import files
import utils
import plot_utils



def funct(reg_dataframe_path : str, cooler_path : str, resolution : int, output_path : str = "./output/plots/noise_plot.png"):

    #TODO: Either make multiprocess or reduce the amount of cooler.fetch calls that are made.

    reg_dataframe = files.load_dataframe(reg_dataframe_path)

    cooler_object = cooler.Cooler(cooler_path)

    temp_dict = {"chrom":[],"prom_name":[],"enh_name":[],"bin1_start":[],"bin1_end":[],"bin2_start":[],"bin2_end":[],"count":[],"noise":[]}

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

        fetch1_start = bin1_start - resolution
        fetch1_end = bin1_end + resolution
        fetch2_start = bin2_start - resolution
        fetch2_end = bin2_end + resolution

        string1 = f"{chrom}:{fetch1_start}-{fetch1_end}"
        string2 = f"{chrom}:{fetch2_start}-{fetch2_end}"

        # * This line is accountable for 99% of used time
        submatrix = cooler_object.matrix(balance=False).fetch(string1,string2)

        noise = np.sum(submatrix)

        # temp_dict = {"prom_name":[],"enh_name":[],"bin1_start":[],"bin1_end":[],"bin2_start":[],"bin2_end":[],"count":[],"noise":[]}
        temp_dict["chrom"].append(chrom)
        temp_dict["prom_name"].append(prom_name)
        temp_dict["enh_name"].append(enh_name)
        temp_dict["bin1_start"].append(bin1_start)
        temp_dict["bin1_end"].append(bin1_end)
        temp_dict["bin2_start"].append(bin2_start)
        temp_dict["bin2_end"].append(bin2_end)
        temp_dict["count"].append(count)
        temp_dict["noise"].append(noise)

    noise_dataframe = pd.DataFrame(temp_dict)

    files.save_dataframe(noise_dataframe,output_path)

    return noise_dataframe


def plot_noise(noise_dataframe : pd.DataFrame, lesser_res_dataframe : pd.DataFrame, output_path = False):
    ""

    noise_dataframe = pd.DataFrame(noise_dataframe).sort_values(by=["chrom","bin1_start","bin2_start"])

    # Remove duplicate pixels. This will remove some pairs, but we only care about specific pixels right now. 
    duplicates = noise_dataframe.loc[:,["bin1_start","bin1_end","bin2_start","bin2_end"]].duplicated(keep="first",)
    
    noise_dataframe = noise_dataframe[~duplicates].sort_values(by=["chrom","bin1_start","bin2_start"])

    chrom_array = np.array(noise_dataframe["chrom"])
    count_array = np.array(noise_dataframe["count"])
    noise_array = np.array(noise_dataframe["noise"])
    #count_and_noise_array = count_array + noise_array
    noise_count_ratio_array = noise_array / count_array

    import matplotlib.pyplot as plt #type:ignore


    # create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12,10))

    # create a set of unique labels
    unique_labels = set(chrom_array)
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    viridis = plt.cm.get_cmap('gist_rainbow', len(unique_labels))
    color_choices = viridis(np.linspace(0, 1, 24))
    #print(newcolors)
    #exit()
    #color_choices = ['#FFC300', '#FF5733', '#C70039', '#900C3F', '#581845', '#0074D9', '#7FDBFF', '#39CCCC', '#3D9970', '#2ECC40', '#FF851B', '#FF4136', '#B10DC9', '#85144b', '#111111', '#AAAAAA', '#DDDDDD', '#FFFFFF', '#F012BE', '#01FF70', '#FFDC00', '#7FDBFF', '#B10DC9', '#F012BE']
    color_dict = {label: color_choices[i] for i, label in enumerate(unique_labels)}
    #color_dict = {label: plt.cm.gist_rainbow(i) for i, label in enumerate(unique_labels)}
    colors = [color_dict[label] for label in chrom_array]
    label_int = np.unique(chrom_array, return_inverse=True)[1]
    #exit()
    # * Plot 1

    for label in unique_labels:
        #indice = chrom_array == label
        indices = [i for i, x in enumerate(chrom_array) if x == label]
        #print(indices)
        # create a new list with the values from array2 at the indices found above
        uq_count_array = [count_array[i] for i in indices]
        uq_noise_array = [noise_array[i] for i in indices]
        scatter = axs[0, 0].scatter(uq_count_array, uq_noise_array, alpha = 0.4, label=label ,color=color_dict[label], s = 20)


    #scatter = axs[0, 0].scatter(count_array, noise_array, c=colors, s = 20)
    axs[0, 0].set_title('Count vs noise per pair')
    axs[0, 0].set_xlabel('Center frequency')
    axs[0, 0].set_ylabel('Surrounding frequency')
    plt.rcParams['legend.fontsize'] = 'small'
    #plt.legend(loc='upper left', bbox_to_anchor=(0, 1), ncol=3)
    axs[0, 0].legend(loc='upper center', columnspacing=0.2, bbox_to_anchor=(0.22, 1), ncol=3, prop={'size': 6})


    # * Plot 2
    box = axs[0, 1].boxplot([count_array, noise_array], labels=["Count","Noise"])
    axs[0, 1].set_yscale('log')
    #axs[0, 1].set_title()



    #* Plot 3
    # Label by chromsome as int range(0,24)
    label_int = np.unique(chrom_array, return_inverse=True)[1]
    count_sums = np.bincount(label_int, weights=count_array)
    noise_sums = np.bincount(label_int, weights=noise_array)
    noise_averages = count_sums/8



    width = 0.25
    x = np.arange(len(count_sums))

    rects1 = axs[1, 0].bar(x - width, count_sums, width, label='Count')
    rects2 = axs[1, 0].bar(x, noise_sums, width, label='Noise')
    rects3 = axs[1, 0].bar(x + width, noise_averages, width, label='Mean Noise', color = 'purple')

    # add labels for the x and y axes
    axs[1, 0].set_xticks(x)
    axs[1, 0].set_xticklabels(np.unique(chrom_array),rotation=45, ha='right')

    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlabel('Label')
    axs[1, 0].set_ylabel('Sum')

    plt.rcParams['legend.fontsize'] = 'small'
    #plt.legend(loc='upper left', bbox_to_anchor=(0, 1), ncol=3)
    #axs[0, 0].legend(loc='upper center', columnspacing=0.2, bbox_to_anchor=(0.22, 1), ncol=3, prop={'size': 6})

    axs[1, 0].legend(prop={'size': 6})
    # * Plot 4
    ratio_sums = np.bincount(label_int, weights=noise_count_ratio_array)
    count = np.bincount(label_int)
    ratio_average = ratio_sums / count

    width = 0.5

    rects3 = axs[1, 1].bar(x, ratio_average, width, label='Noise per count', color="green")
    axs[1, 1].set_xticks(x)
    axs[1, 1].set_xticklabels(np.unique(chrom_array),rotation=45, ha='right')
    axs[1, 1].set_yscale('log')
    axs[1, 1].legend()

    # Show the plot
    plt.show()
    if output_path: 
        plt.savefig(output_path)
        plot_utils.show_plot(output_path)


    



def generate_bias(cooler_file_path, output_path):


    # Load the cooler file
    c = cooler.Cooler(cooler_file_path)

    # Perform ICE normalization
    bias_vector, converged = cooler.ice.iterative_correction(c,cis_only=True, ignore_diags=2)

    # Save the bias values to a file
    fragments = c.bins()[:]
    fragments['bias'] = bias_vector

    print(bias_vector)
    print(fragments)

    fragments.to_csv(output_path, sep='\t', columns=['chrom', 'start', 'end', 'bias'], index=False)
