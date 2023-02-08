import pytest  # type:ignore
import sys
import matplotlib.scale as scl  #type ignore
sys.path.insert(0, '..')
from Functions import *
from Integration_Functions import *
from Cooler_Plotter import *

test_folder = "./testInput/"

def test_collect_counts_for_chr3():
    bed_file = test_folder + "2022-11-14_Reg-Elements_hESC-H1_Grepped_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    modle_file_path = test_folder + "2022-11-14_Modle-Out-hg38-GRCh38-HM12878_barriers.bed.cool"
    mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"

    chrom_name = "chr3"
    start = "9_700_000"
    end = "9_900_000"

    dataframe = create_list_of_chrom(modle_file_path, bed_file, chrom_name, start, end, mcool_file_path)

    collect_counts_and_plot_two_matrixes(dataframe, start, end, outpath="BIOS5410_FUN_9_700_000_to_9_900_000")
    pass


def test_logScale():
    matrix1 = np.array([[1,5,10],
                        [3, 9, 22],
                        [2,1,4]])


    out_filepath = "test_logScale.png"

    plotter = CoolerPlotter()


    plotter.simple_logScale(matrix1, out_filepath = "test_logScale.png")
        
    #print("Log", scl.LogScale(matrix1))


def mini_test():
    bed_file = test_folder + "2022-11-14_Reg-Elements_hESC-H1_Grepped_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"
    chrom_name = "chr19"
    start = "400_000"
    end = "1_000_000"
    dataframe = create_df_of_cool(mcool_file_path + "::/resolutions/1000", bed_file, chrom_name, 1000, start, end)

    print(dataframe)

    outpath = chrom_name + "-" + start + "-" + end + "-mcool_1000"
    collect_count_and_plot_one_matrix(dataframe, start, end, outpath)

def test_collect_counts_and_plot_two_matrixes_density10():
    bed_file = test_folder + "2022-11-14_Reg-Elements_hESC-H1_Grepped_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    modle_file_path = test_folder + "2022-11-20_Modle-Out-hg38-GRCh38-GM12878-RAD21-barriers-tcd-10.cool"
    mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"

    chrom_name = "chr19"
    start = "400_000"
    end = "3_400_000"
    end = "1_000_000"

    dataframe = create_list_of_chrom(modle_file_path, bed_file, chrom_name, start, end, mcool_file_path)

    outname = chrom_name + "-" + start + "-" + end + "-" + "density10"
    collect_counts_and_plot_two_matrixes(dataframe, start, end, outname)


def validate_collected_matrix():
    bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"    

    modle_tcd1_file_path = test_folder + "OUTPUT.cool"

    chrom_name = "chr19"
    start = "400_000"
    #end = "3_400_000"
    end = "1_000_000"


    mcool_file_path = test_folder + "4DNFI9GMP2J8.mcool"
    mcool_1000_file_path = mcool_file_path + "::/resolutions/1000"
    mcool_5000_file_path = mcool_file_path + "::/resolutions/5000"

    modle_tcd1_dataframe = create_df_of_cool(modle_tcd1_file_path, bed_file_new, chrom_name, 5000, start, end)
    mcool_1000_dataframe = create_df_of_cool(mcool_1000_file_path, bed_file_new, chrom_name, 1000, start, end)
    mcool_1000_downressed_dataframe = collect_dataframe_to_match_resolution(mcool_1000_dataframe, target_resolution=5000)
    mcool_5000_dataframe = create_df_of_cool(mcool_5000_file_path, bed_file_new, chrom_name, 5000, start, end)

    plotter = CoolerPlotter()

    modle_tcd1_matrix = plotter.dataframe_to_matrix(modle_tcd1_dataframe, start, end)
    mcool_1000_matrix = plotter.dataframe_to_matrix(mcool_1000_dataframe, start, end)
    mcool_1000_downressed_matrix = plotter.dataframe_to_matrix(mcool_1000_downressed_dataframe, start, end)

    for row in mcool_1000_matrix:
        pass

        
    for row_index, row in enumerate(mcool_1000_matrix):
        if row_index % 5 != 0: continue
        if row_index >= len(mcool_1000_matrix) - 1: continue
        for col_index, column in enumerate(row):
            if col_index % 5 != 0: continue
            if col_index >= len(row) - 1: continue

            sliced = mcool_1000_matrix[row_index:row_index+5, col_index:col_index+5]
            #if on diagonal, delete lower left part

            if col_index == row_index:

                acc = 1
                for i in range(1,5):
                    for y in range(0,acc):
                        sliced[i][y] = 0
                    acc += 1


            sum_of_block = np.sum(sliced)
            small_res_value = mcool_1000_downressed_matrix[row_index // 5][col_index // 5]
            if sum_of_block != small_res_value: 
                print("Incorrect sum")
                print("1000_matrix val:", sum_of_block)
                print("5000_matrix val:", small_res_value)
                print("sliced matrix:\n", sliced)
                print("indexes: ", row_index, col_index)
                print("small_indexes", row_index/5, col_index/5)
                
                exit()
    
    print("Matrix is correct")

def compare_highest_counts():
    bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"    
    modle_tcd1_file_path = test_folder + "2022-11-14_Modle-Out-hg38-GRCh38-HM12878_barriers.bed.cool"
    modle_tcd1_file_path = test_folder + "OUTPUT.cool"

    #mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"
    mcool_file_path = test_folder + "4DNFI9GMP2J8.mcool"
    mcool_1000_file_path = mcool_file_path + "::/resolutions/1000"
    mcool_5000_file_path = mcool_file_path + "::/resolutions/5000"

    chrom_name = "chr19"
    start = "400_000"
    #end = "3_400_000"
    end = "1_000_000"

    modle_tcd1_dataframe = create_df_of_cool(modle_tcd1_file_path, bed_file_new, chrom_name, 5000, start, end)
    mcool_1000_dataframe = create_df_of_cool(mcool_1000_file_path, bed_file_new, chrom_name, 1000, start, end)
    mcool_1000_downressed_dataframe = collect_dataframe_to_match_resolution(mcool_1000_dataframe, target_resolution=5000)
    mcool_5000_dataframe = create_df_of_cool(mcool_5000_file_path, bed_file_new, chrom_name, 5000, start, end)

    plotter = CoolerPlotter()

    modle_tcd1_matrix = plotter.dataframe_to_matrix(modle_tcd1_dataframe, start, end)
    mcool_1000_matrix = plotter.dataframe_to_matrix(mcool_1000_dataframe, start, end)
    mcool_1000_downressed_matrix = plotter.dataframe_to_matrix(mcool_1000_downressed_dataframe, start, end)
    
    modle_tcd1_matrix, mcool_1000_matrix, mcool_1000_downressed_matrix = plotter.equalize_matrices([modle_tcd1_matrix, mcool_1000_matrix, mcool_1000_downressed_matrix])

    modle_max = np.amax(modle_tcd1_matrix)
    mcool_1000_max = np.amax(mcool_1000_downressed_matrix)
    mcool_5000_max = np.amax(mcool_1000_matrix)

    print(modle_max, mcool_1000_max, mcool_5000_max)

def test_create_distance_matrix():
    bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"    
    modle_tcd10_file_path = test_folder + "OUTPUT.cool"

    chrom_name = "chr19"
    start = "400_000"
    #end = "3_400_000"
    end = "1_000_000"

    modle_tcd10_dataframe = create_df_of_cool(modle_tcd10_file_path, bed_file_new, chrom_name, 5000, start, end)

    new_dataframe = create_distance_matrix(modle_tcd10_dataframe)

    print(new_dataframe)

def test_data_based_on_distance():
    bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"    
    modle_tcd10_file_path = test_folder + "OUTPUT.cool"

    chrom_name = "chr19"
    start = "400_000"
    #end = "3_400_000"
    end = "1_000_000"

    modle_tcd10_dataframe = create_df_of_cool(modle_tcd10_file_path, bed_file_new, chrom_name, 5000, start, end)

    new_dataframe = dataframe_based_on_distance(modle_tcd10_dataframe)

    print("New_dataframe:\n ", new_dataframe)


def generate_plots_for_validation_analyzis():
    
    #bed_file = test_folder + "2022-11-14_Reg-Elements_hESC-H1_Grepped_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"    
    modle_tcd1_file_path = test_folder + "2022-11-14_Modle-Out-hg38-GRCh38-HM12878_barriers.bed.cool"
    modle_tcd10_file_path = test_folder + "2022-11-20_Modle-Out-hg38-GRCh38-GM12878-RAD21-barriers-tcd-10.cool"
    modle_tcd1_file_path = test_folder + "OUTPUT.cool"
    modle_tcd10_file_path = test_folder + "OUTPUT.cool"

    #mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"
    mcool_file_path = test_folder + "4DNFI9GMP2J8.mcool"
    mcool_1000_file_path = mcool_file_path + "::/resolutions/1000"
    mcool_5000_file_path = mcool_file_path + "::/resolutions/5000"

    chrom_name = "chr19"
    start = "400_000"
    #end = "3_400_000"
    end = "1_000_000"



    modle_tcd1_dataframe = create_df_of_cool(modle_tcd1_file_path, bed_file_new, chrom_name, 5000, start, end)
    modle_tcd10_dataframe = create_df_of_cool(modle_tcd10_file_path, bed_file_new, chrom_name, 5000, start, end)
    mcool_1000_dataframe = create_df_of_cool(mcool_1000_file_path, bed_file_new, chrom_name, 1000, start, end)
    mcool_1000_downressed_dataframe = collect_dataframe_to_match_resolution(mcool_1000_dataframe, target_resolution=5000)
    mcool_5000_dataframe = create_df_of_cool(mcool_5000_file_path, bed_file_new, chrom_name, 5000, start, end)

    print("\n\n\n\nDataframe RAW: ")
    print(modle_tcd1_dataframe)
    print(modle_tcd10_dataframe)
    print(mcool_1000_dataframe)
    print(mcool_5000_dataframe)


    plotter = CoolerPlotter()

    modle_tcd10_distance_dataframe = create_distance_matrix(modle_tcd10_dataframe)
    plotter.coorelate_matrix(modle_tcd10_distance_dataframe)
    mcool_5000_distance_dataframe = create_distance_matrix(mcool_5000_dataframe)
    plotter.coorelate_matrix(mcool_5000_distance_dataframe)

    exit()

    modle_tcd1_matrix = plotter.dataframe_to_matrix(modle_tcd1_dataframe, start, end)
    modle_tcd10_matrix = plotter.dataframe_to_matrix(modle_tcd10_dataframe, start, end)
    mcool_1000_matrix = plotter.dataframe_to_matrix(mcool_1000_dataframe, start, end)
    mcool_5000_matrix = plotter.dataframe_to_matrix(mcool_5000_dataframe, start, end)
    mcool_1000_downressed_matrix = plotter.dataframe_to_matrix(mcool_1000_downressed_dataframe, start, end)

    #the_shape = modle_tcd1_matrix.flatten().shape
    # distance_matrix_modle = np.empty()
    # for index_row, row in enumerate(modle_tcd1_matrix):
    #     for index_col, col in enumerate(row):
    #         if index_col > index_row: continue
    #         distance = index_col - index_row
    #         distance_index = next(i for i,v in enumerate(distance_matrix_modle) if v.contains(str(distance)))
    #         if distance_index:
    #             distance_matrix_modle[distance_index][1] += col
    #         else:
    #             distance_matrix_modle.append((str(distance),col))

    # print(sum(c for d,c in distance_matrix_modle))

    # the_shape = mcool_1000_matrix.flatten().shape
    # distance_matrix_mcool = np.zeros(the_shape)
    # for index_row, row in enumerate(mcool_1000_matrix):
    #     for index_col, col in enumerate(row):
    #         if index_col > index_row: continue
    #         distance = index_col - index_row
    #         distance_matrix_mcool[distance] += col

    # distance_matrix = [distance_matrix_modle, distance_matrix_mcool]

    # distance_dataframe = pd.DataFrame(distance_matrix)
    # print(distance_dataframe)
    # cormat = distance_dataframe.corr()
    # round(cormat,2)

    exit()

    outname = chrom_name + "-" + start + "-" + end
    outname_modle_tcd1 = outname + "-" + "modle_tcd1.png"
    outname_modle_tcd10 = outname + "-" + "modle_tcd10.png"
    outname_mcool_1000 = outname + "-" + "mcool1000.png"
    outname_mcool_1000_downressed = outname + "-" + "mcool1000_downressed.png"
    outname_mcool_5000 = outname + "-" + "mcool5000.png"

    outname_scatterplot = outname + "-scatterplot.png"
    outname_scatterplot2 = outname + "-scatterplot2.png"

    #print(mcool_1000_dataframe.to_csv("dataframe1.csv"))
    #print(mcool_1000_downressed_dataframe.to_csv("dataframe2.csv"))

    modle_tcd1_matrix, modle_tcd10_matrix, mcool_1000_matrix, mcool_1000_downressed_matrix, mcool_5000_matrix = plotter.equalize_matrices([modle_tcd1_matrix, modle_tcd10_matrix, mcool_1000_matrix, mcool_1000_downressed_matrix, mcool_5000_matrix])


    plotter.scatter_plot(modle_tcd1_matrix, mcool_5000_matrix, outname_scatterplot, open_in_viewer = True)
    plotter.scatter_plot(modle_tcd1_matrix, mcool_1000_downressed_matrix, outname_scatterplot2, open_in_viewer = True)

    plotter.simple_matrix_plot(modle_tcd1_matrix, outname_modle_tcd1, True, axis_start = int(start), axis_end = int(end))
    plotter.simple_matrix_plot(modle_tcd10_matrix, outname_modle_tcd10, True, axis_start = int(start), axis_end = int(end))
    plotter.simple_matrix_plot(mcool_1000_matrix, outname_mcool_1000, True, axis_start = int(start), axis_end = int(end))
    plotter.simple_matrix_plot(mcool_1000_downressed_matrix, outname_mcool_1000_downressed, True, axis_start = int(start), axis_end = int(end))
    plotter.simple_matrix_plot(mcool_5000_matrix, outname_mcool_5000, True, axis_start = int(start), axis_end = int(end))

    #coef_array = plotter.correlation_coeffecient(modle_tcd1_matrix, mcool_5000_matrix)

    coef_array = np.corrcoef(modle_tcd1_matrix.flatten(),mcool_5000_matrix.flatten())
    flattened_modle_tcd1_matrix = modle_tcd1_matrix.flatten()
    flattened_mcool_5000_matrix = mcool_5000_matrix.flatten()


    print(coef_array.shape)
    print("sizes:",len(coef_array.flatten()),len(flattened_modle_tcd1_matrix), len(flattened_mcool_5000_matrix))
    exit()
    for index, ele in enumerate(coef_array.flatten()):
        if ele != 0 and ele != False: 
            print(ele)#, flattened_mcool_5000_matrix[index], flattened_modle_tcd1_matrix[index])
            

    #plotter.simple_matrix_plot(coef_array, "coef_array.png", True, axis_start = int(start), axis_end = int(end))
    exit()


def test_collect_counts_and_plot_two_matrixes():
    bed_file = test_folder + "2022-11-14_Reg-Elements_hESC-H1_Grepped_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ.7group.bed"
    modle_file_path = test_folder + "2022-11-14_Modle-Out-hg38-GRCh38-HM12878_barriers.bed.cool"
    mcool_file_path = test_folder + "2022-11-14_MicroC_hESC-H1_4DNFI9GMP2J8.mcool"

    chrom_name = "chr19"
    start = "400_000"
    end = "3_400_000"
    end = "1_000_000"



    dataframe = create_list_of_chrom(modle_file_path, bed_file, chrom_name, start, end, mcool_file_path)

    print(dataframe)
    exit()
    outname = chrom_name + "-" + start + "-" + end
    collect_counts_and_plot_two_matrixes(dataframe, start, end, outname)


    # dataframe = create_list_of_chrom(modle_file_path, bed_file, chrom_name, start, end, mcool_file_path)

    # outname = chrom_name + "-" + start + "-" + end
    # collect_counts_and_plot_two_matrixes(dataframe, start, end, outname)
    pass

def test_normalize_matrix():
    matrix = np.array([[1,2,3],[3,2,1],[1,1,1]])

    print(matrix)
    print(normalize_matrix(matrix))


def test_equal_matrix():
    matrix = np.array([[1,2,3],[3,2,1],[1,1,1]])

def plot_experiment():
    import matplotlib.pyplot as plt
    import matplotlib
    matrix1 = np.array([[1,5,10],
                        [3, 9, 22],
                        [2,1,4]])
    matrix2 = np.array([[10,2,4],
                        [25,21,25],
                        [3,8,2]])

    axis_start = 300_000
    axis_end = 1_000_000

    matrix_shape = matrix1.shape
    x_ticks_num = matrix_shape[0]
    x_ticks_num = 3
    y_ticks_num = matrix_shape[1]

    axis_total = axis_end - axis_start
    axis_interval = round(axis_total / (x_ticks_num -1))

    x_ticks_array = np.array([])
    for i in range(x_ticks_num):
        tick = axis_total - (axis_interval * i) + axis_start
        print(tick)
        x_ticks_array = np.append(x_ticks_array, tick)
    x_ticks_array = np.flip(x_ticks_array)

    print(x_ticks_array)
    exit()
    ## 0 - 700_000
    axis_total = axis_end - axis_start
    axis_interval = round(axis_total / x_ticks_num)

    x_ticks_array = np.array([])
    for i in range(x_ticks_num + 1):
        tick = axis_total - (axis_interval * i)
        print(tick)
        x_ticks_array = np.append(x_ticks_array, tick)

    print(x_ticks_array)
    exit()
    #plt.xticks([10, 100,200,300,400,500])
    #plt.set_xscale(1,'linear')
    # matplotlib.axes.Axes.set_xscale(1, 'linear')
    #plt.yticks([-1,-2,-3,-4,-5])
    #plt.axis([100,200,300,400])
    #fig, ax0 = plt.subplots()
    #ax0.set_xscale('linear')
    plt.scatter(matrix1,matrix2)
    #max(matrix1)
    plt.xticks(range(0, np.amax(matrix1)), [i for i in range(10, np.amax(matrix1) + 10)])
    plt.savefig("plot_xlab_test.png")
    #plt.arrange()
    plotter = CoolerPlotter()
    plotter.view_plot("plot_xlab_test.png")
    
    exit()
    fig, ax0 = plt.subplots()
    ax0.set_xscale('log')
    formatter0 = EngFormatter(unit ='Hz')
    ax0.xaxis.set_major_formatter(formatter0)
    ax0.plot(xs, ys)
    ax0.set_xlabel('Frequency')
    
    fig.suptitle('matplotlib.axes.Axes.set_xscale() \
    function Example\n', fontweight ="bold")
  
plt.show()

def main():
    pass
    #test_plot_dataframes()

main()