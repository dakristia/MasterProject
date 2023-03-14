import correlation_between_resolutions
import files

def main():
    print("ayo")
#"output/"
    input_folder = "../../input/"
    output_folder = "../../output/"

    #* Modle FIles
    dataframe1_file_name = "cache/chromosome_wide_promoter_enhancer_data_outfile_binsize5000_tcd7.31.cool_H1-hESC.7group.bed.chr1.csv"
    dataframe2_file_name = "chromosome_wide_promoter_enhancer_data_outfile_binsize1000_tcd2.33.cool_H1-hESC.7group.bed.csv"

    #* H1 files

    dataframe1_file_name = f"{input_folder}/dataframes/H1_binsize5000/genome_wide/5000_H1-hESC.7group.bed"
    dataframe2_file_name = f"{input_folder}/dataframes/H1_binsize1000/genome_wide/1000_H1-hESC.7group.bed"


    dataframe1_filepath = output_folder + dataframe1_file_name
    dataframe2_filepath = output_folder + dataframe2_file_name

    dataframe1 = files.load_dataframe(dataframe1_filepath)
    dataframe2 = files.load_dataframe(dataframe2_filepath)
    
    correlation_between_resolutions.box_plots(dataframe1, dataframe2)

    correlation_between_resolutions.correlation_between_resolution(dataframe1,dataframe2)

if __name__ == "__main__":
    main()