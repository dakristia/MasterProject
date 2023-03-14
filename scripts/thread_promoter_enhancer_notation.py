from generate_plots import *
import files

def main():
    test_dataframes()

def test_dataframes(cooler_file_name = "outfile_binsize5000_tcd7.31.cool",
                    bed_file_name = "H1-hESC.7group.bed",
                    chrom_name = "chr19",
                    start : int = False,
                    end : int = False,
                    max_distance : int = 3_000_000):
    
    start = 0
    end = 10_000_000

    cooler_file_path = input_folder + cooler_file_name
    bed_file_path = input_folder + bed_file_name

    ## * Threaded run

    start_time = time.time()

    re_interaction_dataframe : pd.DataFrame = pelg.note_promoter_enhancer_interactions_threaded(
        cooler_file_path = cooler_file_path,
        bed_file_path=bed_file_path,
        chrom_name = chrom_name,
        start = start,
        end = end,
        max_distance = max_distance)
    
    logfile_name_threaded_dataframes = "logfile_threaded_dataframes.log"

    save_dataframe(re_interaction_dataframe, "threaded_dataframe_test1.csv")
    end_time = time.time() - start_time
    print()
    files.append_text_to_file(f"Total time for threaded noting: {end_time} seconds",logfile_name_threaded_dataframes)

    start_time = time.time()
    re_interaction_dataframe_old : pd.DataFrame = pelg.filter_cooler_to_regelement_dataframe(
        cool_file_path=cooler_file_path,
        bed_file_path=bed_file_path,
        chrom_name=chrom_name,
        start=start,
        end=end,
        resolution=5000
    )
    re_interaction_dataframe_old
    save_dataframe(re_interaction_dataframe_old, "unthreaded_dataframe_test1.csv")
    end_time = time.time() - start_time
    files.append_text_to_file(f"Total time for NON-threaded noting: {end_time} seconds",logfile_name_threaded_dataframes)

    files.append_text_to_file(f"Length of threaded dataframe: {len(re_interaction_dataframe)}",logfile_name_threaded_dataframes)
    files.append_text_to_file(f"Length of NON-threaded dataframe: {len(re_interaction_dataframe_old)}",logfile_name_threaded_dataframes)

    print(re_interaction_dataframe)
    print(re_interaction_dataframe_old)    
    re_interaction_dataframe_old



    duplicates = pd.concat([re_interaction_dataframe,re_interaction_dataframe_old]).duplicated(keep="last")

    files.append_text_to_file(f"Duplicates present in dataframes: {len(duplicates)}. \n{duplicates}",logfile_name_threaded_dataframes)


if __name__ == "__main__":
    main()