from generate_plots import *




def test_dataframes(cooler_file_name = "outfile_binsize1000_tcd5.cool",
                    bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed",
                    chrom_name = "chr19",
                    start : int = 0,
                    end : int = 10_000_000,
                    max_distance : int = 3_000_000):
    
    cooler_file_path = input_folder + cooler_file_name
    bed_file_path = input_folder + bed_file_name

    ## * Threaded run

    start_time = time.time()
    promoter_enhancer_dataframe = pelg.extract_pls_els_from_bed(bed_file_path)
    
    re_interaction_dataframe : pd.DataFrame = pelg.note_promoter_enhancer_interactions_threaded(
        cooler_file_path = cooler_file_path,
        promoter_enhancer_dataframe = promoter_enhancer_dataframe,
        chrom_name = chrom_name,
        start = start,
        end = end,
        max_distance = max_distance)
    
    logfile_name_threaded_dataframes = "logfile_threaded_dataframes.log"

    save_dataframe(re_interaction_dataframe, "threaded_dataframe_test1.csv")
    end_time = time.time() - start_time
    print()
    append_text_to_file(f"Total time for threaded noting: {end_time}",logfile_name_threaded_dataframes)

    start_time = time.time()
    re_interaction_dataframe_old : pd.DataFrame = pelg.create_regElem_df_of_cool(
        cool_file_path=cooler_file_path,
        bed_file_path=bed_file_path,
        chrom_name=chrom_name,
        start=start,
        end=end,
        resolution=5000
    )
    save_dataframe(re_interaction_dataframe_old, "unthreaded_dataframe_test1.csv")
    end_time = time.time() - start_time
    append_text_to_file(f"Total time for NON-threaded noting: {end_time}",logfile_name_threaded_dataframes)