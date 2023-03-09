from generate_plots import *

def calculate_statistics_of_one_dataset(chrom_name = "chr19",
                            bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                            cooler_file_name = "outfile_binsize1000_tcd5.cool",
                            start = False,
                            end = False,
                            csv_target_file_path=False):
    
    start_time = time.time()
    print("Starting function 'calculate_statistics_of_one_dataset'...")
    #coolerplotter = plotter.CoolerPlotter()

    bed_file_path = input_folder + bed_file_name
    cooler_file_path = input_folder + cooler_file_name

    cooler_matrix_numpyfilename = cooler_file_name + '.npy'
    cooler_distancematrix_numpyfilename = cooler_file_name + '.distancematrix.npy'
    
    # # Acquire resolution
    # c = cooler.Cooler(cooler_file_path)
    # resolution = int(c.info["bin-size"])

    returnmatrix = load_or_create_matrix(re_matrix_filename = cooler_matrix_numpyfilename, 
                    re_distance_matrix_filename=cooler_distancematrix_numpyfilename,
                    chrom_name="chr19",
                    scale = 5000,
                    cooler_file_path=cooler_file_path,
                    bed_file_path=bed_file_path,
                    )

    re_matrix, re_distance_matrix_expanded = returnmatrix[0], returnmatrix[1]

    # * Start calculating statistics

    ## * Calculate bins with above-zero count

    print("Calculating above-zero counts")

    total_counts = 0

    number_above_zero = 0
    number_above_ten = 0
    number_at_zero = 0

    row_max, col_max = re_matrix.shape
    diagonal_bp_width_modle = 3_000_000
    diagonal_bin_width_modle = (diagonal_bp_width_modle / 5000)# + 1

    for row in range(row_max):
        for col in range(col_max):
            if col > row: break # We only care about one side of the matrix, so don't count anything above diagonal
            if col - row > diagonal_bin_width_modle: 
                break;   # MoDLE only counts up to 3Mbp distance. 
            
            cell_count = re_matrix[row][col]
            if cell_count > 10: number_above_ten += 1
            if cell_count > 0: number_above_zero += 1
            else: number_at_zero += 1
            total_counts += cell_count


    logfile_modle_valid_data = "logfile_modle_vs_valid.txt"
    print()
    append_text_to_file(("Statistics for: " + cooler_file_name), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 0: " + str(number_above_zero)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 10: " + str(number_above_ten)), logfile_modle_valid_data)
    append_text_to_file(("Mean of matrix: " + str(np.mean(re_matrix))), logfile_modle_valid_data)
    append_text_to_file(("Mean of matrix above 0: " + str(total_counts / number_above_zero)), logfile_modle_valid_data)

    end_time = time.time() - start_time

    #columns = ["chrom","resolution","pe_total_count","pe_mean_count","pe_std_count","pe_median_count",
                "pe_min_count","pe_max_count","pe_total_bins","chr_total_bins","pe_bin_frequency","nonpe_bin_frequency"]

    dictionary = {"chrom":np.array([]),"resolution":np.array([]),"pe_total_count":np.array([]),"pe_mean_count":np.array([]),"pe_std_count":np.array([]),
                "pe_median_count":np.array([]),"pe_min_count":np.array([]), "pe_max_count":np.array([]),"pe_total_bins":np.array([]),"chr_total_bins":np.array([]),
                "pe_bin_frequency":np.array([]),"nonpe_bin_frequency":np.array([])}
    


    # dictionary = {'chrom':[chrom_name],'resolution':[resolution],'total_bins_in_chrom':[number_of_bins],
    #             'total_bins_with_pls_and_els':[total_bins_with_pls_and_els],
    #             'min_count':[np.min(counter_np)],'max_count':[np.max(counter_np)],'average_count':[average_pls_els_per_bin],
    #             'median_count':[np.median(counter_np)],'standarddeviation':[np.std(counter_np)],'total_count':[np.sum(counter_np)],
    #             'list_of_counts':[counter_np],'list_of_indexes':[coordinate_np]}

def calculate_statistics(chrom_name = "chr19",
                            bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                            modle_5000_file_name = "outfile_binsize5000_tcd10.cool", 
                            #modle_1000_file_name = "outfile_binsize1000_tcd10.cool",
                            modle_1000_file_name = "outfile_binsize1000_tcd5.cool",
                            valid_25000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/25000",
                            valid_5000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/5000",
                            valid_1000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/1000",
                            start = "400_000",
                            end = "1_000_000"):
    
    start = "1_000_000"
    end = "1_600_000"

    start_time = time.time()
    print("Starting function 'calculate_statistics'...")
    coolerplotter = plotter.CoolerPlotter()

    bed_file_path = input_folder + bed_file_name
    modle_5000_file_path = input_folder + modle_5000_file_name
    modle_1000_file_path = input_folder + modle_1000_file_name
    valid_25000_file_path = input_folder + valid_25000_file_name
    valid_5000_file_path = input_folder + valid_5000_file_name
    valid_1000_file_path = input_folder + valid_1000_file_name


    ## * Acquire cooler-objects and their respective 2D selectors
    cooler_object_time_start = time.time()
    print("Acquiring cooler objects")
    modle_5000_cool_object = pelg.File_Functions.read_mcool_file(modle_5000_file_path)
    modle_5000_cool_2D_selector = modle_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    modle_1000_cool_object = pelg.File_Functions.read_mcool_file(modle_1000_file_path)
    modle_1000_cool_2D_selector = modle_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    valid_25000_cool_object = pelg.File_Functions.read_mcool_file(valid_25000_file_path)
    valid_25000_cool_2D_selector = valid_25000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    valid_5000_cool_object = pelg.File_Functions.read_mcool_file(valid_5000_file_path)
    valid_5000_cool_2D_selector = valid_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    valid_1000_cool_object = pelg.File_Functions.read_mcool_file(valid_1000_file_path)
    valid_1000_cool_2D_selector = valid_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    ## * Define filenames
    
    modle_5000_numpyfilename = f'{modle_5000_file_name}.{chrom_name}.{start}.{end}'
    modle_1000_numpyfilename = f'{modle_1000_file_name}.{chrom_name}.{start}.{end}'
    valid_5000_numpyfilename = f'{valid_5000_file_name}.{chrom_name}.{start}.{end}'.replace('::/resolutions/','')
    valid_1000_numpyfilename = f'{valid_1000_file_name}.{chrom_name}.{start}.{end}'.replace('::/resolutions/','')
    modle_5000_matrix_re_numpyfilename = f'{modle_5000_numpyfilename}.matrix.re.npy'
    modle_1000_matrix_re_numpyfilename = f'{modle_1000_numpyfilename}.matrix.re.npy'
    valid_5000_matrix_re_numpyfilename = f'{valid_5000_numpyfilename}.matrix.re.npy'
    valid_1000_matrix_re_numpyfilename = f'{valid_1000_numpyfilename}.matrix.re.npy'

    modle_5000_distance_matrix_re_numpyfilename = f'{modle_5000_numpyfilename}.distancematrix.re.npy'
    modle_1000_distance_matrix_re_numpyfilename = f'{modle_1000_numpyfilename}.distancematrix.re.npy'
    valid_5000_distance_matrix_re_numpyfilename = f'{valid_5000_numpyfilename}.distancematrix.re.npy'
    valid_1000_distance_matrix_re_numpyfilename = f'{valid_1000_numpyfilename}.distancematrix.re.npy'

    # * Defining variables to make them available in function scope
    re_modle_5000_dataframe = None
    re_modle_5000_matrix = None
    re_modle_1000_dataframe = None
    re_modle_1000_dataframe_scaled_5000 = None
    re_modle_1000_matrix_scaled_5000 = None
    re_valid_5000_dataframe = None
    re_valid_5000_matrix = None
    re_valid_1000_dataframe = None
    re_valid_1000_dataframe_scaled_5000 = None
    re_valid_1000_matrix_scaled_5000 = None

    re_dataframe_start_time = time.time()
    print("Fetching regulatory element dataframe in:",re_dataframe_start_time)

    ## * Regulatory elements MoDLE data res 5000

    modle_5000_matrices = load_or_create_matrix(re_matrix_filename=modle_5000_matrix_re_numpyfilename,
                                                re_distance_matrix_filename=modle_5000_distance_matrix_re_numpyfilename,
                                                cooler_file_path=modle_5000_file_path,
                                                bed_file_path=bed_file_path,
                                                chrom_name=chrom_name)



    re_modle_5000_matrix = modle_5000_matrices[0]
    re_modle_5000_distance_matrix_expanded = modle_5000_matrices[1]

    modle_1000_matrices = load_or_create_matrix(re_matrix_filename=modle_1000_matrix_re_numpyfilename,
                                                re_distance_matrix_filename=modle_1000_distance_matrix_re_numpyfilename,
                                                cooler_file_path=modle_1000_file_path,
                                                bed_file_path=bed_file_path,
                                                chrom_name=chrom_name)



    re_modle_1000_matrix = modle_1000_matrices[0]
    re_modle_1000_distance_matrix_expanded = modle_1000_matrices[1]

    valid_5000_matrices = load_or_create_matrix(re_matrix_filename=valid_5000_matrix_re_numpyfilename,
                                                re_distance_matrix_filename=valid_5000_distance_matrix_re_numpyfilename,
                                                cooler_file_path=valid_5000_file_path,
                                                bed_file_path=bed_file_path,
                                                chrom_name=chrom_name)



    re_valid_5000_matrix = valid_5000_matrices[0]
    re_valid_5000_distance_matrix_expanded = valid_5000_matrices[1]

    valid_1000_matrices = load_or_create_matrix(re_matrix_filename=valid_1000_matrix_re_numpyfilename,
                                                re_distance_matrix_filename=valid_1000_distance_matrix_re_numpyfilename,
                                                cooler_file_path=valid_1000_file_path,
                                                bed_file_path=bed_file_path,
                                                chrom_name=chrom_name)



    re_dataframe_end_time = time.time() - re_dataframe_start_time
    print("Finished fetching regulatory element dataframes in:",re_dataframe_end_time)

    first_column = [x[0] for x in re_modle_5000_distance_matrix_expanded]



    ## * Percentage misses
    # For every element in 1000 res (scaled to 5000) array: 
    #   If any pixel at position a,b is bigger than 0, and the pixel in the same position in re_modle_5000_matrix
    #   is 0, increase number_of_missed by 1. 

    total_counts_5000 = 0
    total_counts_1000 = 0

    number_above_zero_5000 = 0
    number_above_zero_1000 = 0
    number_above_ten_5000 = 0
    number_above_ten_1000 = 0
    number_at_zero_5000 = 0
    number_at_zero_1000 = 0

    row_max, col_max = re_modle_1000_matrix_scaled_5000.shape
    diagonal_bp_width_modle = 3_000_000
    diagonal_bin_width_modle = (diagonal_bp_width_modle / 5000)# + 1

    for row in range(row_max):
        for col in range(col_max):
            if col > row: break # We only care about one side of the matrix, so don't count anything above diagonal
            if col - row > diagonal_bin_width_modle: 
                if cell_count_re_5000 > 0: print(re_modle_5000_matrix[row][col]); exit();
                break;   # MoDLE only counts up to 3Mbp distance. 
            
            cell_count_re_5000 = re_modle_5000_matrix[row][col]
            if cell_count_re_5000 > 10: number_above_ten_5000 += 1
            if cell_count_re_5000 > 0: number_above_zero_5000 += 1
            else: number_at_zero_5000 += 1
            total_counts_5000 += cell_count_re_5000

            cell_count_re_1000 = re_modle_1000_matrix_scaled_5000[row][col]
            if cell_count_re_1000 > 10: number_above_ten_1000 += 1
            if cell_count_re_1000 > 0: number_above_zero_1000 += 1
            else: number_at_zero_1000 += 1
            total_counts_1000 += cell_count_re_1000



    difference_above_zero = number_above_zero_5000  - number_above_zero_1000
    difference_ten_zero = number_above_ten_5000  - number_above_ten_1000
    #percentage_false_positive = difference_above_zero / number_above_zero_5000 * 100

    logfile_modle_valid_data = "logfile_modle_vs_valid.txt"
    print()
    append_text_to_file("MoDLE statistics.", logfile_modle_valid_data)
    append_text_to_file(("Numbers above 0 at 5000bp res: " + str(number_above_zero_5000)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 0 at 1000bp res: " + str(number_above_zero_1000)), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 5000:" + str(np.mean(re_modle_5000_matrix))), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 1000:" + str(np.mean(re_modle_1000_matrix_scaled_5000))), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 5000 above 0:" + str(total_counts_5000 / number_above_zero_5000)), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 1000 above 0:" + str(total_counts_1000 / number_above_zero_1000)), logfile_modle_valid_data)
 
    append_text_to_file(("Difference above 0: " + str(difference_above_zero)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 10 at 5000bp res: " + str(number_above_ten_5000)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 10 at 1000bp res: " + str(number_above_ten_1000)), logfile_modle_valid_data)    
    append_text_to_file(("Difference above 10: " + str(difference_ten_zero)), logfile_modle_valid_data)    
    #write_to_file(("Percent difference: " + str(percentage_false_positive)), logfile_modle_valid_data)

    number_above_zero_5000 = 0
    number_above_zero_1000 = 0
    number_above_ten_5000 = 0
    number_above_ten_1000 = 0

    number_at_zero_5000 = 0
    number_at_zero_1000 = 0

    row_max, col_max = re_valid_1000_matrix_scaled_5000.shape
    for row in range(row_max):
        for col in range(col_max):
            if col > row: break # We only care about one side of the matrix, so don't count anything above diagonal
            if col - row > diagonal_bin_width_modle:  break # MoDLE only counts up to 3Mbp distance. 
            
            cell_count_re_5000 = re_valid_5000_matrix[row][col]
            if cell_count_re_5000 > 10: number_above_ten_5000 += 1
            if cell_count_re_5000 > 0: number_above_zero_5000 += 1
            else: number_at_zero_5000 += 1
            total_counts_5000 += cell_count_re_5000

            cell_count_re_1000 = re_valid_1000_matrix_scaled_5000[row][col]
            if cell_count_re_1000 > 10: number_above_ten_1000 += 1
            if cell_count_re_1000 > 0: number_above_zero_1000 += 1
            else: number_at_zero_1000 += 1
            total_counts_1000 += cell_count_re_1000


    difference_above_zero = number_above_zero_5000  - number_above_zero_1000
    difference_ten_zero = number_above_ten_5000  - number_above_ten_1000
    print()
    append_text_to_file("Valid statistics.", logfile_modle_valid_data)
    append_text_to_file(("Numbers above 0 at 5000bp res: " + str(number_above_zero_5000)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 0 at 1000bp res: " + str(number_above_zero_1000)), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 5000:" + str(np.mean(re_valid_5000_matrix))), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 1000:" + str(np.mean(re_valid_1000_matrix_scaled_5000))), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 5000 above 0:" + str(total_counts_5000 / number_above_zero_5000)), logfile_modle_valid_data)
    append_text_to_file(("Mean of array 1000 above 0:" + str(total_counts_1000 / number_above_zero_1000)), logfile_modle_valid_data)
    
    append_text_to_file(("Difference above 0: " + str(difference_above_zero)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 10 at 5000bp res: " + str(number_above_ten_5000)), logfile_modle_valid_data)
    append_text_to_file(("Numbers above 10 at 1000bp res: " + str(number_above_ten_1000)), logfile_modle_valid_data)        

    append_text_to_file(("Difference above 10: " + str(difference_ten_zero)), logfile_modle_valid_data)

    #write_to_file(("Percent difference: " + str(percentage_false_positive)), logfile_modle_valid_data)

    ## * Correlation coeffecient

    re_modle_5000_matrix = re_modle_5000_matrix / re_modle_5000_matrix.sum()
    re_valid_5000_matrix = re_valid_5000_matrix / re_valid_5000_matrix.sum()
    re_modle_1000_matrix_scaled_5000 = re_modle_1000_matrix_scaled_5000 / re_modle_1000_matrix_scaled_5000.sum()
    re_valid_1000_matrix_scaled_5000 = re_valid_1000_matrix_scaled_5000 / re_valid_1000_matrix_scaled_5000.sum()
    
    matrices = [re_modle_5000_matrix , re_valid_5000_matrix, re_modle_1000_matrix_scaled_5000, re_valid_1000_matrix_scaled_5000]

    for row in range(len(matrices)):
        for j in range(row+1, len(matrices)):
            corr_coeff = np.corrcoef(matrices[row].flatten(), matrices[j].flatten())[0,1]
            print(f"Correlation coefficient between matrix {row+1} and matrix {j+1}: {corr_coeff}")

    f,ax = plt.subplots(figsize=(7,6))
    plt.title("Correlation Matrix")


    matrices = [x.flatten() for x in matrices]
    corr_matrix = np.corrcoef(matrices, rowvar=False)
    print(corr_matrix)
    plt.imshow(corr_matrix, cmap='YlOrRd', interpolation='nearest')
    plt.colorbar()

    plt.savefig("corr_image.png")
    coolerplotter.view_plot("corr_image.png")

    end_time = time.time() - start_time
    print("Function 'calculate_statistics' finished runtime in:", end_time, "s")
