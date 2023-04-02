import files
import pandas as pd
import numpy as np
import concurrent.futures
import plot_utils
from matplotlib import pyplot as plt # type:ignore

def _get_correlating_counts(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame):
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

            #  
            handled_bins_dataframe = pd.concat((handled_bins_dataframe,bins_in_range_dataframe),ignore_index=True)

            df1_count_index = dataframe1.columns.get_loc('modle_count') + 1
            df1_count = row[df1_count_index]

            df2_unique_bins = bins_in_range_dataframe.groupby(['bin1_start', 'bin1_end', 'bin2_start', 'bin2_end']).agg('first').reset_index()

            # if len(df2_unique_bins) > 1:
            #     print(bins_in_range_dataframe)
            #     print(df2_unique_bins)
            #     exit()
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

            if df2_coarsen_count > df1_count and df1_count > 2:
            
                print(df1_count)
                print(df2_coarsen_count)
                print(df2_counts)
                print(row)
                print(bins_in_range_dataframe)
                exit()

            df1_original_counts = np.append(df1_original_counts,df1_count)
            df2_coarsen_counts = np.append(df2_coarsen_counts,df2_coarsen_count)


    return df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsen_counts

def plot_correlation_between_resolution(dataframe1 : pd.DataFrame, dataframe2 : pd.DataFrame, label1: str, label2 : str, output_path : str, 
                                    res1 : int = False, res2 : int = False, logScale : bool = False, scatter_plot : bool = True, box_plot : bool = True,
                                    workers : int = 5):
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
            print(f)
            print("hjosada")
            returned_df1, returned_df2, returned_original_df1, returned_coarsened_df2 = f.result()
            
            df1_flatten_counts = np.append(df1_flatten_counts,returned_df1)
            df2_flatten_counts = np.append(df2_flatten_counts,returned_df2) 
            df1_original_counts = np.append(df1_original_counts,returned_original_df1)
            df2_coarsened_counts = np.append(df2_coarsened_counts,returned_coarsened_df2) 

    if scatter_plot and box_plot: 
        fig = plt.figure(figsize=(12, 7))
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        ax1 = plt.subplot2grid((2, 3), (0, 0))
        ax2 = plt.subplot2grid((2, 3), (0, 1))
        ax3 = plt.subplot2grid((2, 3), (1, 0))
        ax4 = plt.subplot2grid((2, 3), (1, 1))
        ax5 = plt.subplot2grid((2, 3), (0, 2), colspan=2, rowspan=2)


    else: fig, ax1 = plt.subplots()

    if logScale: plt.yscale('log'); plt.xscale('log') 

    if scatter_plot:
        
        ax1.scatter(df1_flatten_counts, df2_flatten_counts, c=df1_flatten_counts, cmap = 'hsv')
        ax1.set_xlabel(f'raw counts \n{res1} bp resolution')
        ax1.set_ylabel(f'raw counts \n{res2} bp resolution')
        ax1.set_title(f'A')

        ## * With correlation line
        ax2.scatter(df1_flatten_counts, df2_flatten_counts, c=df1_flatten_counts, cmap = 'hsv')
        ax2.set_xlabel(f'raw counts \n{res1} bp resolution')
        ax2.set_ylabel(f'raw counts \n{res2} bp resolution')
        ax2.set_title(f'B')

        min_val = np.min([df1_flatten_counts.min(), df2_flatten_counts.min()])
        max_val = np.max([df1_flatten_counts.max(), df2_flatten_counts.max()])

        perfect_correlation_x = np.linspace(min_val,max_val)
        perfect_correlation_y = perfect_correlation_x
        ax2.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")
        
        ## * Plot original df1 with coarsened df2
        ax3.scatter(df1_original_counts, df2_coarsened_counts, c=df1_original_counts, cmap = 'hsv')
        ax3.set_xlabel(f'raw counts \n{res1} bp resolution')
        ax3.set_ylabel(f'raw counts \n{res2} bp resolution')
        ax3.set_title(f'C')

        ## * With correlation line 
        ax4.scatter(df1_original_counts, df2_coarsened_counts, c=df1_original_counts, cmap = 'hsv')
        ax4.set_xlabel(f'raw counts \n{res1} bp resolution')
        ax4.set_ylabel(f'raw counts \n{res2} bp resolution')
        ax4.set_title(f'D')

        min_val = np.min([df1_original_counts.min(), df2_coarsened_counts.min()])
        max_val = np.max([df1_original_counts.max(), df2_coarsened_counts.max()])

        perfect_correlation_x = np.linspace(min_val,max_val)
        perfect_correlation_y = perfect_correlation_x
        ax4.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")
        
    if box_plot:
        ax5.boxplot([df1_original_counts,df2_flatten_counts], labels=[label1,label2])
        ax5.set_xlabel('Data Source')
        ax5.set_ylabel('Raw Counts')
        ax5.set_title(f'D')

    plt.savefig(output_path)
    plt.show()


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