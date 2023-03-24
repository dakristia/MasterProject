import files
import Plotter as plotter


def compare_tcd_plots(chrom_name = "chr19",
                    bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                    tcd1_filename = "outfile_binsize1000.cool" , 
                    tcd5_filename = "outfile_binsize1000_tcd5.cool",
                    tcd10_filename = "outfile_binsize1000_tcd10.cool",
                    valid_filename = "4DNFI9GMP2J8.mcool::/resolutions/1000",
                    start = "400_000", 
                    end = "1_000_000"):

    input_folder = "../../input/dataframes/distance_dataframes/tcd_comparison_distance_dataframes/"

    bed_file_path = input_folder + bed_file_name

    tcd1_filepath = input_folder + tcd1_filename
    tcd1_numpyfilename = f'{tcd1_filename}.{chrom_name}.{start}.{end}'
    # tcd1_matrix_raw_numpyfilename = f'{tcd1_numpyfilename}.matrix.raw.npy'
    # tcd1_matrix_re_numpyfilename = f'{tcd1_numpyfilename}.matrix.re.npy'
    # tcd1_dataframe_raw_numpyfilename = f'{tcd1_numpyfilename}.dataframe.raw.npy'
    # tcd1_dataframe_re_numpyfilename = f'{tcd1_numpyfilename}.dataframe.re.npy'
    tcd1_distance_dataframe_raw_numpyfilepath = f'{input_folder}{tcd1_numpyfilename}.distancedataframe.raw.npy'
    tcd1_distance_dataframe_re_numpyfilepath = f'{input_folder}{tcd1_numpyfilename}.distancedataframe.re.npy'

    tcd5_filepath = input_folder + tcd5_filename
    tcd5_numpyfilename = f'{tcd5_filename}.{chrom_name}.{start}.{end}'
    # tcd5_matrix_raw_numpyfilename = f'{tcd5_numpyfilename}.matrix.raw.npy'
    # tcd5_matrix_re_numpyfilename = f'{tcd5_numpyfilename}.matrix.re.npy'
    # tcd5_dataframe_raw_numpyfilename = f'{tcd5_numpyfilename}.dataframe.raw.npy'
    # tcd5_dataframe_re_numpyfilename = f'{tcd5_numpyfilename}.dataframe.re.npy'
    tcd5_distance_dataframe_raw_numpyfilepath = f'{input_folder}{tcd5_numpyfilename}.distancedataframe.raw.npy'
    tcd5_distance_dataframe_re_numpyfilepath = f'{input_folder}{tcd5_numpyfilename}.distancedataframe.re.npy'

    tcd10_filepath = input_folder + tcd10_filename
    tcd10_numpyfilename = f'{tcd10_filename}.{chrom_name}.{start}.{end}'
    # tcd10_matrix_raw_numpyfilename = f'{tcd10_numpyfilename}.matrix.raw.npy'
    # tcd10_matrix_re_numpyfilename = f'{tcd10_numpyfilename}.matrix.re.npy'
    # tcd10_dataframe_raw_numpyfilename = f'{tcd10_numpyfilename}.dataframe.raw.npy'
    # tcd10_dataframe_re_numpyfilename = f'{tcd10_numpyfilename}.dataframe.re.npy'
    tcd10_distance_dataframe_raw_numpyfilepath = f'{input_folder}{tcd10_numpyfilename}.distancedataframe.raw.npy'
    tcd10_distance_dataframe_re_numpyfilepath = f'{input_folder}{tcd10_numpyfilename}.distancedataframe.re.npy'

    valid_filepath = input_folder + valid_filename
    valid_numpyfilename = f'{valid_filename}.{chrom_name}.{start}.{end}'.replace("::/resolutions/","")
    # valid_matrix_raw_numpyfilename = f'{valid_numpyfilename}.matrix.raw.npy'
    # valid_matrix_re_numpyfilename = f'{valid_numpyfilename}.matrix.re.npy'
    # valid_dataframe_raw_numpyfilename = f'{valid_numpyfilename}.dataframe.raw.npy'
    # valid_dataframe_re_numpyfilename = f'{valid_numpyfilename}.dataframe.re.npy'
    valid_distance_dataframe_raw_numpyfilepath = f'{input_folder}{valid_numpyfilename}.distancedataframe.raw.npy'
    valid_distance_dataframe_re_numpyfilepath = f'{input_folder}{valid_numpyfilename}.distancedataframe.re.npy'

    coolerplotter = plotter.CoolerPlotter()

    # return_values = load_or_create_matrix(re_matrix_filename=tcd1_matrix_re_numpyfilename,
    #                         cache_dataframe=True,
    #                         bed_file_path=bed_file_path,
    #                         re_dataframe_filepath=tcd1_dataframe_re_numpyfilename,
    #                         chrom_name=chrom_name,
    #                         start=start,
    #                         end=end,
    #                         cooler_file_path=tcd1_filepath)
    # print(return_values)
    # matrix = return_values[0]


    # tcd1_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd1_distance_dataframe_raw_numpyfilename,
    #                                     re_distance_dataframe_filepath= tcd1_distance_dataframe_re_numpyfilename,
    #                                     cooler_file_path=tcd1_filepath,
    #                                     bed_file_path=bed_file_path,
    #                                     chrom_name=chrom_name,
    #                                     start = start,
    #                                     end = end,)
    


    # tcd5_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd5_distance_dataframe_raw_numpyfilename,
    #                                     re_distance_dataframe_filepath= tcd5_distance_dataframe_re_numpyfilename,
    #                                     cooler_file_path=tcd5_filepath,
    #                                     bed_file_path=bed_file_path,
    #                                     chrom_name=chrom_name,
    #                                     start = start,
    #                                     end = end,)
    
    # tcd10_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd10_distance_dataframe_raw_numpyfilename,
    #                                     re_distance_dataframe_filepath= tcd10_distance_dataframe_re_numpyfilename,
    #                                     cooler_file_path=tcd10_filepath,
    #                                     bed_file_path=bed_file_path,
    #                                     chrom_name=chrom_name,
    #                                     start = start,
    #                                     end = end,)

    # valid_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= valid_distance_dataframe_raw_numpyfilename,
    #                                 re_distance_dataframe_filepath= valid_distance_dataframe_re_numpyfilename,
    #                                 cooler_file_path=valid_filepath,
    #                                 bed_file_path=bed_file_path,
    #                                 chrom_name=chrom_name,
    #                                 start = start,
    #                                 end = end,)

    tcd1_distance_raw_dataframe = files.load_dataframe(tcd1_distance_dataframe_raw_numpyfilepath) 
    tcd1_distance_re_dataframe = files.load_dataframe(tcd1_distance_dataframe_re_numpyfilepath)

    tcd5_distance_raw_dataframe = files.load_dataframe(tcd5_distance_dataframe_raw_numpyfilepath) 
    tcd5_distance_re_dataframe = files.load_dataframe(tcd5_distance_dataframe_re_numpyfilepath)

    tcd10_distance_raw_dataframe = files.load_dataframe(tcd10_distance_dataframe_raw_numpyfilepath) 
    tcd10_distance_re_dataframe = files.load_dataframe(tcd10_distance_dataframe_re_numpyfilepath)

    valid_distance_raw_dataframe = files.load_dataframe(valid_distance_dataframe_raw_numpyfilepath) 
    valid_distance_re_dataframe = files.load_dataframe(valid_distance_dataframe_re_numpyfilepath)


    ## * Plot average distance between bins

    coolerplotter.average_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd10",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd5",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd1",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="H1 1000bp",
                                        show = True, plotname="Average counts based distance between bins")
    
    ## * Plot average distance only between bins containing regulatory elements
    coolerplotter.average_distance_plot(dataframe = tcd10_distance_re_dataframe,line_lable="MoDLE 1000bp:tcd10",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(dataframe = tcd5_distance_re_dataframe,line_lable="MoDLE 1000bp:tcd5",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = tcd1_distance_re_dataframe,line_lable="MoDLE 1000bp:tcd1",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = valid_distance_re_dataframe,line_lable="H1 1000bp",
                                        show = True, plotname="Average counts based distance between bins with regulatory elements")


    ## * Plot total counts of all distances
    coolerplotter.totalcounts_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd10",
                                        newplot=True,show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd5",
                                        show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd1",
                                        show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="H1 1000bp",
                                        show = True, plotname="Total counts based distance between bins with regulatory elements")


    ## * Plot total counts of all distances
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd10",
                                        newplot=True,show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd5",
                                        show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd1",
                                        show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="H1 1000bp",
                                        show = True, plotname="Logarithmic standard deviation for regulatory elements with distance")



    coolerplotter.standard_deviation_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd10",
                                        newplot=True,show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd5",
                                        show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000bp:tcd1",
                                        show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="H1 1000bp",
                                        show = True, plotname="Standard deviation for regulatory elements with distance", logY=False)    

    

#     tcd1_matrices = load_or_create_matrix(raw_matrix_filename=tcd1_matrix_raw_numpyfilename,
#                                             re_matrix_filename=tcd1_matrix_re_numpyfilename,
#                                             cooler_file_path=tcd1_filepath,
#                                             chrom_name=chrom_name,
#                                             start=start,
#                                             end=end,
#                                             bed_file_path=bed_file_path,
#                                             raw_distance_dataframe_filepath=tcd1_distance_raw_dataframe,
#                                             re_distance_dataframe_filepath=tcd1_distance_re_dataframe
#                                         )
    
#     tcd1_raw_matrix = tcd1_matrices[0]
#     tcd1_re_matrix = tcd1_matrices[1]

#     tcd5_matrices = load_or_create_matrix(raw_matrix_filename=tcd5_matrix_raw_numpyfilename,
#                                             re_matrix_filename=tcd5_matrix_re_numpyfilename,
#                                             cooler_file_path=tcd5_filepath,
#                                             chrom_name=chrom_name,
#                                             start=start,
#                                             end=end,
#                                             bed_file_path=bed_file_path,
#                                             raw_distance_dataframe_filepath=tcd5_distance_raw_dataframe,
#                                             re_distance_dataframe_filepath=tcd5_distance_re_dataframe
#                                         )
    
#     tcd5_raw_matrix = tcd5_matrices[0]
#     tcd5_re_matrix = tcd5_matrices[1]

#     tcd10_matrices = load_or_create_matrix(raw_matrix_filename=tcd10_matrix_raw_numpyfilename,
#                                             re_matrix_filename=tcd10_matrix_re_numpyfilename,
#                                             cooler_file_path=tcd10_filepath,
#                                             chrom_name=chrom_name,
#                                             start=start,
#                                             end=end,
#                                             bed_file_path=bed_file_path,
#                                             raw_distance_dataframe_filepath=tcd10_distance_raw_dataframe,
#                                             re_distance_dataframe_filepath=tcd10_distance_re_dataframe
#                                         )
    
#     tcd10_raw_matrix = tcd10_matrices[0]
#     tcd10_re_matrix = tcd10_matrices[1]

#     valid_matrices = load_or_create_matrix(raw_matrix_filename=valid_matrix_raw_numpyfilename,
#                                             re_matrix_filename=valid_matrix_re_numpyfilename,
#                                             cooler_file_path=valid_filepath,
#                                             chrom_name=chrom_name,
#                                             start=start,
#                                             end=end,
#                                             bed_file_path=bed_file_path,
#                                             raw_distance_dataframe_filepath=valid_distance_raw_dataframe,
#                                             re_distance_dataframe_filepath=valid_distance_re_dataframe
#                                         )
    
#     valid_raw_matrix = valid_matrices[0]
#     valid_re_matrix = valid_matrices[1]

#     tcd1_standarddeviation = tcd1_distance_raw_dataframe["standarddeviation"].to_numpy()
#     tcd5_standarddeviation = tcd5_distance_raw_dataframe["standarddeviation"].to_numpy()
#     tcd10_standarddeviation = tcd10_distance_raw_dataframe["standarddeviation"].to_numpy()

# #    coolerplotter.scatter_plot(tcd1_standarddeviation,tcd5_standarddeviation,"tcd1","tcd5","Sandarddeviation correlation","standarddeviation_scatterplot.png",open_in_viewer=True)

#     matrices = [tcd1_raw_matrix,tcd5_raw_matrix,tcd10_raw_matrix]
#     for row1 in range(len(matrices)):
#         for row2 in range(row1+1, len(matrices)):
#             corr_coeff = np.corrcoef(matrices[row1].flatten(), matrices[row2].flatten())[0,1]
#             print(f"Correlation coefficient between raw matrix {row1+1} and raw matrix {row2+1}: {corr_coeff}")

#     for row in range(len(matrices)):
#         corr_coeff_valid = np.corrcoef(matrices[row].flatten(), valid_raw_matrix.flatten())[0,1]
#         print(f"Correlation coefficient between raw matrix {row+1} and raw valid matrix: {corr_coeff_valid}")


#     matrices = [tcd1_re_matrix,tcd5_re_matrix,tcd10_re_matrix]
#     for row1 in range(len(matrices)):
#         for row2 in range(row1+1, len(matrices)):
#             corr_coeff = np.corrcoef(matrices[row1].flatten(), matrices[row2].flatten())[0,1]
#             print(f"Correlation coefficient between re matrix {row1+1} and re matrix {row2+1}: {corr_coeff}")


#     for row in range(len(matrices)):
#         corr_coeff_valid = np.corrcoef(matrices[row].flatten(), valid_re_matrix.flatten())[0,1]
#         print(f"Correlation coefficient between re matrix {row+1} and re valid matrix: {corr_coeff_valid}")

