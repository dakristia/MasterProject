import files
import pandas as pd
import numpy as np
import os
import concurrent.futures
import plot_utils
from matplotlib import pyplot as plt # type:ignore
    

def _get_correlating_counts(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame):
    """
    Flattens dataframe1 and dataframe2. For each count in dataframe1, find correlating counts in dataframe2. 
    For each of these counts, copy the single count from dataframe1. Resulting flatten counts will be longer than originals.

    Also create a list of dataframe1's original counts, and create a matching array of dataframe2 by "coarsing" the 
    counts to a lower resolution. E.g. collect all counts for bins that are within the range of bins from dataframe1.


    Args:
        dataframe1 (pd.DataFrame): _description_
        dataframe2 (pd.DataFrame): _description_

    Returns:
        _type_: _description_
    """

    df1_flatten_counts = np.array([])
    df2_flatten_counts = np.array([])
    df1_original_counts = np.array([])
    df2_coarsen_counts = np.array([])

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

            # * Find all bins that correlate 
            bins_in_range_dataframe = dataframe2.query(f'bin1_start >= {bin1_start} and bin1_end <= {bin1_end}') 
            bins_in_range_dataframe = bins_in_range_dataframe.query(f'bin2_start >= {bin2_start} and bin2_end <= {bin2_end}' )  

            
            handled_bins_dataframe = pd.concat((handled_bins_dataframe,bins_in_range_dataframe),ignore_index=True)

            df1_count_index = dataframe1.columns.get_loc('modle_count') + 1
            df1_count = row[df1_count_index]

            df2_unique_bins = bins_in_range_dataframe.groupby(['bin1_start', 'bin1_end', 'bin2_start', 'bin2_end']).agg('first').reset_index()

            # * Get all counts relating to the bin
            df2_counts = np.array(df2_unique_bins['modle_count'])


            # * To calculate correlation, we need equal size arrays.
            # * Duplicate count in lower res array to match the number of found counts in the higher res array
            df1_counts = [df1_count for c in df2_counts]



            if len(df2_counts) == 0:
                df1_counts = np.array([df1_count])
                df2_counts = np.array([0])

            df1_flatten_counts = np.append(df1_flatten_counts,df1_counts)
            df2_flatten_counts = np.append(df2_flatten_counts,df2_counts)

            # * Summarise all counts in relevant high res region
            # * "Coarses" the counts to a lowere resolution

            df2_coarsen_count = np.sum(df2_counts)

            df1_original_counts = np.append(df1_original_counts,df1_count)
            df2_coarsen_counts = np.append(df2_coarsen_counts,df2_coarsen_count)

    return df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsen_counts

def flatten_and_coarsen_data(dataframe1 : pd.DataFrame, dataframe2 : pd.DataFrame,
                                    res1 : int = False, res2 : int = False, workers : int = 5, 
                                    cache = False, cache_path : str = "./cache/flattened_counts_cache.npy"):
    
    # * Extract name and extension from original path
    file_name = os.path.splitext(cache_path)[0]
    file_extension = os.path.splitext(cache_path)[1]

    # * Give filename to each array cache
    df1_flatten_filepath = file_name + "df1flatten" + file_extension
    df2_flatten_filepath = file_name + "df2flatten" + file_extension
    df1_original_filepath = file_name + "df1original" + file_extension
    df2_coarsened_filepath = file_name + "df2coarsen" + file_extension

    # * Get cached data if it exists
    if cache: 
        all_file_found = True
        try:
            df1_flatten_counts = np.load(df1_flatten_filepath)
        except FileNotFoundError:
            print(f"No file found at {df1_flatten_filepath}. Proceeding as normal.")
            all_file_found = False
        try:
            df2_flatten_counts = np.load(df2_flatten_filepath)
        except FileNotFoundError:
            print(f"No file found at {df2_flatten_filepath}. Proceeding as normal.")
            all_file_found = False
        try:
            df1_original_counts = np.load(df1_original_filepath)
        except FileNotFoundError:
            print(f"No file found at {df1_original_filepath}. Proceeding as normal.")
            all_file_found = False
        try:
            df2_coarsened_counts = np.load(df2_coarsened_filepath)
        except FileNotFoundError:
            print(f"No file found at {df2_coarsened_filepath}. Proceeding as normal.")
            all_file_found = False

        if all_file_found: 
            print(f"All cache files found. Returning.")
            return df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsened_counts
        else:
            print(f"At least one cache file missing.")
        #TODO: CHange so function only creates arrays that are missing. RIght now, function will create all arrays even if only one is missing.

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
    df1_original_counts = np.array([])
    df2_coarsened_counts = np.array([])

    with concurrent.futures.ProcessPoolExecutor(max_workers = workers) as executor:
        futures = []
        for chrom in chrom_names:
            chrom_dataframe1 = dataframe1[dataframe1['chrom'] == chrom]
            chrom_dataframe2 = dataframe2[dataframe2['chrom'] == chrom]

            futures.append(executor.submit(_get_correlating_counts, chrom_dataframe1, chrom_dataframe2))
        
        concurrent.futures.wait(futures)

        for f in futures:
            returned_df1, returned_df2, returned_original_df1, returned_coarsened_df2 = f.result()
            
            df1_flatten_counts = np.append(df1_flatten_counts,returned_df1)
            df2_flatten_counts = np.append(df2_flatten_counts,returned_df2)
            df1_original_counts = np.append(df1_original_counts,returned_original_df1)
            df2_coarsened_counts = np.append(df2_coarsened_counts,returned_coarsened_df2)

    if cache:
        np.save(df1_flatten_filepath,df1_flatten_counts)
        np.save(df2_flatten_filepath,df2_flatten_counts)
        np.save(df1_original_filepath,df1_original_counts)
        np.save(df2_coarsened_filepath,df2_coarsened_counts)

    return df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsened_counts

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