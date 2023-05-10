import files

import predict_chromosome_data
import calculate_general_statistics_of_cool_data
import PromoterEnhancerListGenerator as pelg


import os
import pandas as pd
import time



def main():
    #save_matrix()
    #create_dir_for_path("./output/lmao/haha")


    start_time = time.time()

    chrom_name = "chr19"

    #cooler_file_name = "outfile_binsize5000_tcd7.31.cool"
    #cooler_file_name = "outfile_binsize1000_tcd2.33.cool"
    #cooler_file_name = "H1-hESC.mcool:://resolutions/5000"
    cooler_file_name = "H1-hESC.mcool:://resolutions/1000"

    cooler_file_path = "../../input/" + cooler_file_name
    bed_file_name = "H1-hESC.7group.bed"
    bed_file_path = "../../input/" + bed_file_name
    chrom_sizes_file_path = "../../input/" + "hg38.chrom.sizes"
    output_path = f"../../output/chromosome_wide_promoter_enhancer_data_{cooler_file_name}_{bed_file_name}"
    output_path_statistics = f"../../output/statistics_chromosome_wide_promoter_enhancer_data_{cooler_file_name}_{bed_file_name}"


    chromosome_wide_dataframe = pelg.note_promoter_enhancer_interactions_genomewide(cooler_file_path=cooler_file_path,
                                                        bed_file_path=bed_file_path,
                                                        output_path=output_path,
                                                        cache_dataframe=True,
                                                        )


    end_time = time.time() - start_time
    print(f"Finished noting interactions in {end_time} seconds.")

    #note_promoter_enhancer_interactions_chromwide

    # if os.path.isfile("../output/test_dataframe_DELETEME.csv"):
    #     interaction_dataframe = files.load_dataframe("../output/test_dataframe_DELETEME.csv")
    # else:
    #     interaction_dataframe = pelg.Dataframe_Functions.note_promoter_enhancer_interactions_multiprocess(cooler_file_path,bed_file_path,chrom_name)
    # files.save_dataframe(interaction_dataframe, "../output/test_dataframe_DELETEME.csv")

    start_time = time.time()

    calculated_statistics_dataframe = calculate_general_statistics_of_cool_data.calculate_statistics_for_dataframe(chromosome_wide_dataframe)

    end_time = time.time() - start_time
    print(f"Finished calculating statistics in {end_time} seconds.")
    
    files.save_dataframe(calculated_statistics_dataframe, output_path_statistics)
    exit()

    chrom_sizes_dataframe = pd.read_csv(chrom_sizes_file_path,delim_whitespace=True,header=None,names=["chrom","size"])
    chrom_size = int(chrom_sizes_dataframe[chrom_sizes_dataframe['chrom'] == chrom_name]['size'])
    resolution = 5000

    

    predict_chromosome_data.calc_reg_pairs_multiprocess("../../input/" + "H1-hESC.7group.bed", chrom_name, chrom_size, 5000, 10, 20, f"./output/predicted_regulatoryinteraction_{chrom_name}_{resolution}")

    
    
    pass

def test_ssa_note_promoter_enhancer_interactions_multiprocess():
    cooler_file_name = "outfile_binsize1000_tcd2.33.cool"

    cooler_file_path = "../../input/" + cooler_file_name
    bed_file_name = "H1-hESC.7group.bed"
    bed_file_path = "../../input/" + bed_file_name

# 1,chr1,1013929,1014264,778570,778919,EH38E2776808,EH38E2776539,1013000,1014000,778000,779000,3,-1
# 2,chr1,1013929,1014264,778570,778919,EH38E2776808,EH38E2776539,1014000,1015000,778000,779000,3,-1
    returned_df = pelg.note_promoter_enhancer_interactions_multiprocess(cooler_file_path,bed_file_path,"chr1",0,1_200_000,workers=1)
    #returned_df = pelg._note(, bed_file_path, "chr1", 0, 1_100_000, 0)
    print(returned_df)
    exit()


if __name__ == "__main__":
    main()