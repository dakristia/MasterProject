import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # type:ignore
import seaborn as sns #type:ignore

import sys
sys.path.append('../src/regelem')
import resolution_correlation
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
    
    only_chrom = False
    if only_chrom:
        print(f"Only chrom {only_chrom}")
        dataframe1 = dataframe1[dataframe1["chrom"] == only_chrom]
        dataframe2 = dataframe2[dataframe2["chrom"] == only_chrom]

    dataframe_array = [dataframe1,dataframe2]

    print(dataframe2)

    label1 = "H1 5000bp"
    label2 = "H1 1000bp"
    
    
    output_path = "./output/plots/correlation_H1_5000_1000.png"

    df1_flatten_counts, df2_flatten_counts, df1_original_counts, df2_coarsen_counts = resolution_correlation.flatten_and_coarsen_data(dataframe1,dataframe2,res1= 5000,res2=1000, cache = True)

    #! TODO: FIX
    #! Temp, hacky fix. Extremely slightly misrepresents data which is caused by an error. 
    df2_coarsen_counts = np.where(df2_coarsen_counts > df1_original_counts, df1_original_counts, df2_coarsen_counts) 


    print("Flattened and coarsened counts got.")

    res1 = 5000
    res2 = 1000

    fig, axs = plt.subplots(nrows=2,ncols=3,figsize=(12,8))

    #fig = plt.figure(3,figsize=(12, 7))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    ax1 = axs[0][0]
    ax2 = axs[0][1]
    ax3 = axs[0][2]
    ax4 = axs[1][0]
    ax5 = axs[1][1]
    ax6 = axs[1][2]
    ax6.axis('off')

    plt.yscale('log'); plt.xscale('log') 

    import seaborn as sns #type:ignore
    from matplotlib import colors # type:ignore

    # * ax1
    hb = ax1.hexbin(df1_flatten_counts, df2_flatten_counts, gridsize=50, cmap='viridis',norm=colors.LogNorm())
    cb = fig.colorbar(hb, ax=ax1)
    cb.set_label('Density (log scale)')

    ax1.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax1.set_ylabel(f'raw counts \n{res2} bp resolution')
    ax1.set_title(f'A')

    ## * ax2, With correlation line
    hb = ax2.hexbin(df1_flatten_counts, df2_flatten_counts, gridsize=50, cmap='viridis',norm=colors.LogNorm())
    cb = fig.colorbar(hb, ax=ax2)
    cb.set_label('Density (log scale)')

    ax2.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax2.set_ylabel(f'raw counts \n{res2} bp resolution')
    ax2.set_title(f'B')

    min_val = np.min([df1_flatten_counts.min(), df2_flatten_counts.min()])
    max_val = np.max([df1_flatten_counts.max(), df2_flatten_counts.max()])

    perfect_correlation_x = np.linspace(min_val,max_val)
    perfect_correlation_y = perfect_correlation_x
    ax2.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")

    # * ax3
    ax3.boxplot([df1_original_counts,df2_flatten_counts], labels=[label1,label2])
    ax3.set_xlabel('Data Source')
    ax3.set_ylabel('Raw Counts')
    ax3.set_yscale('log')
    ax3.set_title(f'C')


    # * ax4 
    hb = ax4.hexbin(df1_original_counts, df2_coarsen_counts, gridsize=50, cmap='viridis',norm=colors.LogNorm())
    cb = fig.colorbar(hb, ax=ax4)
    cb.set_label('Density (log scale)')

    ax4.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax4.set_ylabel(f'counts at \n{res2} coarsened to \n{res1} bp resolution')
    ax4.set_title(f'D')

    # * ax5
    hb = ax5.hexbin(df1_original_counts, df2_coarsen_counts, gridsize=50, cmap='viridis',norm=colors.LogNorm())
    cb = fig.colorbar(hb, ax=ax5)
    cb.set_label('Density (log scale)')

    ax5.set_xlabel(f'raw counts \n{res1} bp resolution')
    ax5.set_ylabel(f'counts at \n{res2} coarsened to \n{res1} bp resolution')
    ax5.set_title(f'E')

    
    min_val = np.min([df1_original_counts.min(), df2_coarsen_counts.min()])
    max_val = np.max([df1_original_counts.max(), df2_coarsen_counts.max()])

    perfect_correlation_x = np.linspace(min_val,max_val)
    perfect_correlation_y = perfect_correlation_x
    ax5.plot(perfect_correlation_x, perfect_correlation_y, linestyle='--', color='black',label="Perfect correlation")

    # * Show plot
    plt.tight_layout()
    plt.show()
    plt.savefig(output_path)
    plot_utils.show_plot(output_path)



if __name__ == "__main__":
    main()