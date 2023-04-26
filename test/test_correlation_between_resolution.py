import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # type:ignore
import seaborn as sns #type:ignore

import sys
sys.path.append('../src/regelem')
import correlation_between_resolutions
import files
import plot_utils

def main():
    print("ayo")

    input_folder = "../input/"
    output_folder = "../output/"

    # * Modle FIles
    dataframe1_file_name = "cache/chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.chr1.csv"
    dataframe2_file_name = "chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv"

    # * H1 files genome wide

    dataframe1_file_name = f"{input_folder}dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.csv"
    dataframe2_file_name = f"{input_folder}dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.csv"
    
    # # * Test Files
    #dataframe1_file_name = f"{input_folder}dataframes/H1_binsize5000/5000_H1-hESC.7group.bed.chrY.csv"
    #dataframe2_file_name = f"{input_folder}dataframes/H1_binsize1000/1000_H1-hESC.7group.bed.chrY.csv"

    dataframe1_filepath = dataframe1_file_name
    dataframe2_filepath = dataframe2_file_name

    dataframe1 = files.load_dataframe(dataframe1_filepath)
    dataframe2 = files.load_dataframe(dataframe2_filepath)
    
    dataframe_array = [dataframe1,dataframe2]

    print(dataframe2)

    label1 = "H1 5000bp"
    label2 = "H1 1000bp"
    
    
    output_path = "./output/plots/correlation_H1_5000_1000.png"

    df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsen_counts = correlation_between_resolutions.plot_correlation_between_resolution(dataframe1,dataframe2,res1= 5000,res2=1000)

    res1 = 5000
    res2 = 1000


    fig, axs = plt.subplots(nrows=1,ncols=3,figsize=(12,3))

    #fig = plt.figure(3,figsize=(12, 7))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    ax1 = axs[0]#plt.subplot2grid((3), (0))
    ax2 = axs[1]#plt.subplot2grid((3), (1))
    ax5 = axs[2]#plt.subplot2grid((3), (2), colspan=2)

    plt.yscale('log'); plt.xscale('log') 

    import seaborn as sns #type:ignore
    from scipy.stats import gaussian_kde

    #sns.kdeplot(df1_flatten_counts, df2_flatten_counts, cmap="Blues", shade=True, shade_lowest=False)
    # * Claculate the density of the points
    vstacked_flatten_counts = np.vstack([df1_flatten_counts,df2_flatten_counts])
    gaussian = gaussian_kde(vstacked_flatten_counts)(vstacked_flatten_counts)

    # * Sort the points by density, so that the densest points are plotted last
    idx = gaussian.argsort()
    x, y, gaussian = df1_flatten_counts[idx], df2_flatten_counts[idx], gaussian[idx]

    flattened_dataframe = pd.DataFrame({"df1":df1_flatten_counts,"df2":df2_flatten_counts})

    sns.scatterplot(
        data=flattened_dataframe,
        x = "df1",
        y = "df2",
        c=gaussian,
        cmap="viridis",
        ax=ax1,
)

    #ax1.scatter(x, y, c=gaussian, s=50)

    #ax1.scatter(df1_flatten_counts, df2_flatten_counts)
    ax1.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax1.set_ylabel(f'raw counts \n{res2} bp resolution')
    ax1.set_title(f'A')

    ## * With correlation line
    ax2.scatter(x, y, c=gaussian, s=50)
    #ax2.scatter(df1_flatten_counts, df2_flatten_counts)
    ax2.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax2.set_ylabel(f'raw counts \n{res2} bp resolution')
    ax2.set_title(f'B')

    min_val = np.min([df1_flatten_counts.min(), df2_flatten_counts.min()])
    max_val = np.max([df1_flatten_counts.max(), df2_flatten_counts.max()])

    perfect_correlation_x = np.linspace(min_val,max_val)
    perfect_correlation_y = perfect_correlation_x
    ax2.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")
    
    # ## * Plot original df1 with coarsened df2
    # # * Claculate the density of the points
    # vstacked_flatten_counts = np.vstack([df1_original_counts,df2_coarsened_counts])
    # gaussian = gaussian_kde(vstacked_flatten_counts)(vstacked_flatten_counts)

    # # * Sort the points by density, so that the densest points are plotted last
    # idx = gaussian.argsort()
    # x, y, gaussian = df1_original_counts[idx], df2_coarsened_counts[idx], gaussian[idx]

    # ax3.scatter(x, y, c=gaussian, s=50)
    # #ax3.scatter(df1_original_counts, df2_coarsened_counts)
    # ax3.set_xlabel(f'raw counts \n{res1} bp resolution')
    # ax3.set_ylabel(f'raw counts \n{res2} bp resolution')
    # ax3.set_title(f'C')

    # ## * With correlation line
    # ax4.scatter(x, y, c=gaussian, s=50)
    # #ax4.scatter(df1_original_counts, df2_coarsened_counts)
    # ax4.set_xlabel(f'raw counts \n{res1} bp resolution')
    # ax4.set_ylabel(f'raw counts \n{res2} bp resolution')
    # ax4.set_title(f'D')

    # min_val = np.min([df1_original_counts.min(), df2_coarsened_counts.min()])
    # max_val = np.max([df1_original_counts.max(), df2_coarsened_counts.max()])

    # perfect_correlation_x = np.linspace(min_val,max_val)
    # perfect_correlation_y = perfect_correlation_x
    # ax4.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")
    
    
    ax5.boxplot([df1_original_counts,df2_flatten_counts], labels=[label1,label2])
    ax5.set_xlabel('Data Source')
    ax5.set_ylabel('Raw Counts')
    ax5.set_title(f'D')

    plt.show()
    plt.savefig(output_path)
    plot_utils.show_plot(output_path)

if __name__ == "__main__":
    main()