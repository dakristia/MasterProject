
import sys
sys.path.append('../src/regelem')
import read_and_plot_predicted_chromosome_data

def main():
    #read_and_plot_predicted_chromosome_data.plot_something()

    read_and_plot_predicted_chromosome_data.plot_average_regulatory_interaction()
    #ile_path = "../../output/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.csv"
    #read_and_plot_predicted_chromosome_data.reformat_statistics_file_and_get_dict(file_path,5000)

if __name__ == "__main__": 
    main()

    #"/mnt/e/UiO/Master/MasterProject/output/statistics_chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.csv"