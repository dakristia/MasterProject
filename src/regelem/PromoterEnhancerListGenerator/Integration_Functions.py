# import warnings
# import cooler                   #type:ignore
# import matplotlib as mpl        #type:ignore
# #mpl.use('TkAgg')    
# import matplotlib.pyplot as plt #type:ignore
# import numpy as np              #type:ignore
# import pandas as pd             #type:ignore
# from collections import defaultdict

# from .File_Functions import *
# # from Constants import *
# # from File_Functions import *
# # from Functions import *
# from IPython.display import display
# import re
# import time
# import sys
# import os

# #TODO Test this

# absolute_path = os.path.dirname(__file__)
# plotter_path = absolute_path + "/../../Plotter"
# listgen_path = absolute_path + "/../../PromoterEnhancerListGenerator"


# sys.path.insert(0, '../../Plotter')
# #sys.path.insert(0, '../../PromoterEnhancerListGenerator')

# from Cooler_Plotter import * #type:ignore

# def equalize_matrix(matrix : np.array):
#     sumM = np.sum(matrix)

#     equalized_matrix = matrix / sumM

#     return equalized_matrix

# def normalize_matrix(matrix : np.array):
#     highest = np.amax(matrix)

#     normalized_matrix = matrix / highest

#     return normalized_matrix

# def png_duplicate_rename(filepath : str) -> str:
#     """Checks if input feilpath already exists. If yes, renames path. 

#     Args:
#         filepath (str): filepath to check

#     Returns:
#         str: new filepath
#     """
#     if os.path.exists(filepath):
#         filepath = filepath.replace(".png", "") #Temporarity remove file ending
#         filepath = filepath + "_(1).png"

#     dupcounter = 1
#     while os.path.exists(filepath):
#         dupcounter += 1
#         filepath = filepath.replace(".png", "") #Temporarity remove file ending
#         pathlength = len(filepath)
#         replacement = "_(" + str(dupcounter) + ").png"
#         filepath = filepath.replace(filepath[pathlength - 4:], replacement)
    
#     return filepath


# def collect_count_and_plot_one_matrix(dataframe : pd.core.frame.DataFrame, 
#                     start : str = None, 
#                     end : str = None,
#                     outpath : str = None):
#     print("Creating plot for a single dataframe")

#     if outpath == None:
#         outpath = "output/plotMatrix"
#         outpath = outpath + "Modle.png"
#         outpath = png_duplicate_rename(outpath)
#     else:
#         outpath = outpath + "_Modle.png"
#         outpath = png_duplicate_rename(outpath)

#     resolution = abs(dataframe["bin2_start"][0] - dataframe["bin2_end"][0])

#     print(resolution)

#     if start == None: start = min(dataframe["bin2_start"])
#     else: start = int(start)

#     if end == None: end = max(dataframe["bin2_end"])
#     else: end = int(end)

#     print(start,end,resolution)

#     length = int((end-start)/resolution)

#     print(length)

#     modle_matrix = np.zeros(shape=[length,length])
#     validation_matrix = np.zeros(shape=[length,length])

#     print(dataframe['modle_count'])

#     #Only need one of these
#     #Consider making a function for this
#     for row in dataframe.itertuples():
#         row_index = row[0]
#         bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
#         bin1_start = row[bin1_start_index + 1]
#         bin1_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_end")
#         bin1_end = row[bin1_end_index + 1]

#         bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
#         bin2_start = row[bin2_start_index+ 1]
#         bin2_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_end")
#         bin2_end = row[bin2_end_index + 1]

#         matrix_row : int = int((bin1_start - start)/ resolution)
#         matrix_column : int = int((bin2_start - start)/ resolution)

#         modle_count_index = DATAFRAME_COLUMNS_INTERNAL.index("modle_count")
#         modle_count = row[modle_count_index + 1]

#         valid_count_index = DATAFRAME_COLUMNS_INTERNAL.index("valid_count")
#         valid_count = row[valid_count_index + 1]

#         modle_matrix[matrix_row][matrix_column] += modle_count
#         validation_matrix[matrix_row][matrix_column] += valid_count
        
#         if matrix_row != matrix_column:
#             modle_matrix[matrix_column][matrix_row] += modle_count
#             validation_matrix[matrix_column][matrix_row] += valid_count




#     plotter = CoolerPlotter()

#     print(modle_matrix)

#     plotter.simple_matrix_plot(modle_matrix, outpath, True, axis_start = start, axis_end = end)




# def collect_counts_and_plot_two_matrixes(dataframe : pd.core.frame.DataFrame, 
#                     start : str = None, 
#                     end : str = None,
#                     outpath : str = None) -> None:
#     """Reads count columns of dataframe and make two matrices of them. Create interaction plots and produce .png of them"

#     Args:
#         dataframe (pd.core.frame.DataFrame): _description_
#         start (str, optional): _description_. Defaults to None.
#         end (str, optional): _description_. Defaults to None.
#         outpath (str, optional): _description_. Defaults to None.
#     """
#     #TODO Change doc comment to be more specific on columns
    
#     print("Creating plots for dataframes")
    
#     print(dataframe)

#     if outpath == None:
#         outpath = "output/plotMatrix"
#         outpathMatrix1 = outpath + "Modle.png"
#         outpathMatrix1 = png_duplicate_rename(outpathMatrix1)
#         outpathMatrix2 = outpath + "Valid.png"
#         outpathMatrix2 = png_duplicate_rename(outpathMatrix2)
#         outpathCompMatrix = outpath + "Comp.png"
#         outpathCompMatrix = png_duplicate_rename(outpathCompMatrix)
#     else:
#         outpathMatrix1 = outpath + "_Modle.png"
#         outpathMatrix1 = png_duplicate_rename(outpathMatrix1)
#         outpathMatrix2 = outpath + "_Valid.png"
#         outpathMatrix2 = png_duplicate_rename(outpathMatrix2)
#         outpathCompMatrix = outpath + "Comp.png"
#         outpathCompMatrix = png_duplicate_rename(outpathCompMatrix)
        


#     resolution = abs(dataframe["bin2_start"][0] - dataframe["bin2_end"][0])

#     if start == None: start = min(dataframe["bin2_start"])
#     else: start = int(start)

#     if end == None: end = min(dataframe["bin2_end"])
#     else: end = int(end)

#     print(start,end,resolution)

#     length = int((end-start)/resolution)

#     print(length)

#     modle_matrix = np.zeros(shape=[length,length])
#     validation_matrix = np.zeros(shape=[length,length])

#     for row in dataframe.itertuples():
#         row_index = row[0]
#         bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
#         bin1_start = row[bin1_start_index + 1]
#         bin1_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_end")
#         bin1_end = row[bin1_end_index + 1]

#         bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
#         bin2_start = row[bin2_start_index+ 1]
#         bin2_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_end")
#         bin2_end = row[bin2_end_index + 1]

#         matrix_row : int = int((bin1_start - start)/ resolution)
#         matrix_column : int = int((bin2_start - start)/ resolution)

#         modle_count_index = DATAFRAME_COLUMNS_INTERNAL.index("modle_count")
#         modle_count = row[modle_count_index + 1]

#         valid_count_index = DATAFRAME_COLUMNS_INTERNAL.index("valid_count")
#         valid_count = row[valid_count_index + 1]

#         modle_matrix[matrix_row][matrix_column] += modle_count
#         validation_matrix[matrix_row][matrix_column] += valid_count
        
#         if matrix_row != matrix_column:
#             modle_matrix[matrix_column][matrix_row] += modle_count
#             validation_matrix[matrix_column][matrix_row] += valid_count


#     comparisonMatrix = abs(modle_matrix - validation_matrix)

#     #modle_matrix = equalize_matrix(modle_matrix)

#     ratio = np.sum(modle_matrix) / np.sum(validation_matrix)
#     validation_matrix = validation_matrix * ratio

#     print("HERE")
#     print(np.amax(modle_matrix))
#     print(np.sum(modle_matrix))
#     print(np.amax(validation_matrix))
#     print(np.sum(validation_matrix))
#     print("HERE")
    

#     print(modle_matrix.shape)

#     plotter = CoolerPlotter()

#     #modle_matrix = equalize_matrix(modle_matrix)
#     #validation_matrix = equalize_matrix(validation_matrix)

#     plotter.simple_matrix_plot(modle_matrix, outpathMatrix1, True, axis_start = start, axis_end = end)
#     plotter.simple_matrix_plot(validation_matrix, outpathMatrix2, True, axis_start = start, axis_end = end)
#     #plotter.simple_matrix_plot(comparisonMatrix, outpathCompMatrix, True, normalize = False, axis_start = start, axis_end = end)



#     #! Scatter plots

#     xAxis = dataframe["bin1_start"].to_numpy()
#     yAxis = dataframe["bin2_start"].to_numpy()

#     plotter.scatter_plot(xAxis, yAxis)