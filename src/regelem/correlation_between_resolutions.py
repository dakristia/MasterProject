import files
import pandas as pd
import numpy as np
import concurrent.futures
import plot_utils
from matplotlib import pyplot as plt # type:ignore

def _get_correlating_counts(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame):
    df1_flatten_counts = np.array([])
    df2_flatten_counts = np.array([])

    # * DataFrame to collect all rows that we query in the following for-loop. Will later use to check if any were missed. 
    handled_bins_dataframe = pd.DataFrame(columns=dataframe1.columns)

    for row in dataframe1.itertuples():
            # * With itertuples we index on the row instead of by column name
            bin1_start_index = dataframe1.columns.get_loc('bin1_start') + 1
            bin1_end_index = dataframe1.columns.get_loc('bin1_end') + 1
            bin2_start_index = dataframe1.columns.get_loc('bin2_start') + 1
            bin2_end_index = dataframe1.columns.get_loc('bin2_end') + 1

            bin1_start = row[bin1_start_index]
            bin1_end = row[bin1_end_index]
            bin2_start = row[bin2_start_index]
            bin2_end = row[bin2_end_index]

            # * Find all bins that are 
            bins_in_range_dataframe = dataframe2.query(f'bin1_start >= {bin1_start} and bin1_end <= {bin1_end}') 
            bins_in_range_dataframe = bins_in_range_dataframe.query(f'bin2_start >= {bin2_start} and bin2_end <= {bin2_end}' )  

            #  
            handled_bins_dataframe = pd.concat((handled_bins_dataframe,bins_in_range_dataframe),ignore_index=True)

            df1_count_index = dataframe1.columns.get_loc('modle_count') + 1
            df1_count = row[df1_count_index]

            # * Get all counts relating to the bin
            df2_counts = np.array(bins_in_range_dataframe['modle_count'])

            # * To calculate correlation, we need equal size arrays.
            df1_counts = [df1_count for c in df2_counts]

            df1_flatten_counts = np.append(df1_flatten_counts,df1_counts)
            df2_flatten_counts = np.append(df2_flatten_counts,df2_counts)

    

    return df1_flatten_counts, df2_flatten_counts

def correlation_between_resolution(dataframe1 : pd.DataFrame, dataframe2 : pd.DataFrame, output_path : str,
                                    res1 : int = False, res2 : int = False, logScale : bool = False, scatter_plot : bool = True, box_plot : bool = True):
    #TODO: Make this dynamic.
    if res2 > res1:
        raise Exception("res2 higher than res1. Excepted res1 to be higher than res2.")

    if not res1:
        res1 = int(dataframe1['bin1_end'][0] - dataframe1['bin1_start'][0])

    if not res2:
        res2 = int(dataframe2['bin2_end'][0] - dataframe2['bin2_start'][0])


    chrom_names = np.array(pd.unique(dataframe1['chrom']))

    df1_flatten_counts = np.array([])
    df2_flatten_counts = np.array([])


    with concurrent.futures.ProcessPoolExecutor(max_workers = 5) as executor:

        futures = []
        for chrom in chrom_names:
            chrom_dataframe1 = dataframe1[dataframe1['chrom'] == chrom]
            chrom_dataframe2 = dataframe2[dataframe2['chrom'] == chrom]

            futures.append(executor.submit(_get_correlating_counts, chrom_dataframe1, chrom_dataframe2))

        concurrent.futures.wait(futures)

        for f in futures:
            returned_df1, returned_df2 = f.result()
            
            df1_flatten_counts = np.append(df1_flatten_counts,returned_df1)
            df2_flatten_counts = np.append(df2_flatten_counts,returned_df2)

    
    fig, ax = plt.subplots()

    if logScale: plt.yscale('log')
    if scatter_plot:
        plt.scatter(df1_flatten_counts, df2_flatten_counts)
        plt.xlabel(f'Raw contact counts at {res1}bp resolution')
        plt.ylabel(f'Raw contact counts at {res2}bp resolution')
        plt.title(f'Correlation between counts at {res1}bp resolution and {res2}bp resolution.')
        plt.show()
        
        plt.savefig(output_path)

    # if box_plot:
    #     dictionary = {'modle_counts':}
    #     box_plots

#! Todo: Change box_plots to use lists of lists in stead of dataframes
def box_plots(dataframes : list, labels : list):

    # * Fuck numpy arrays
    all_counts = list()
    index = 0

    for df in dataframes:

        counts_array = list(df['modle_count'])
        all_counts.append(counts_array)
        index+=1

    print(len(all_counts))
    print(labels)

    # * Create a boxplot
    fig, ax = plt.subplots(figsize=(8,10))
    ax.boxplot(all_counts,labels=labels)
    # * Set the title and axis labels
    #ax.set_title('')
    ax.set_xlabel('Data Source')
    ax.set_ylabel('Count')
    # * Show the plot
    plt.show()
    plt.savefig("../../output/plots/boxplot_example.png")
    plot_utils.view_plot("../../output/plots/boxplot_example.png")

    ax.set_yscale('log')
    plt.show()
    plt.savefig("../../output/plots/boxplot_example_log.png")
    
    plot_utils.view_plot("../../output/plots/boxplot_example_log.png")

def violin_plot(dataframes : list , labels : list):

    all_counts = list()
    for df in dataframes:

        counts_array = list(df['modle_count'])
        all_counts.append(counts_array)
        

    fig, ax = plt.subplots(
        figsize=(10,5))

    ax.violinplot(all_counts)
    #ax[0].violinplot(all_counts[0],showmeans=True,showextrema=True)
    #ax[1].violinplot(all_counts[1],showmeans=True,showextrema=True)
    # * Set the title and axis labels
    #ax.set_title('')
    #ax.set_xlabel('Data Source')
    #ax.set_ylabel('Raw Counts')
    # * Show the plot
    plt.show()
    plt.savefig("../../output/plots/violinplot_example.png")
    plot_utils.view_plot("../../output/plots/violinplot_example.png")

    ax.set_yscale('log')
    plt.show()
    plt.savefig("../../output/plots/violinplot_example_log.png")
    
    plot_utils.view_plot("../../output/plots/violinplot_example_log.png")