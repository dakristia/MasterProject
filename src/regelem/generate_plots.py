import sys
import os
from os.path import dirname, basename, isfile, join

sys.path.insert(0, '..')
from Constants import *
import numpy as np
import pandas as pd
import PromoterEnhancerListGenerator as pelg
import Plotter as plotter
import matplotlib
import matplotlib.pyplot as plt
import cooler
import time
import math
import psutil
from datetime import datetime
import threading
import concurrent.futures


test_folder = "./testInput/"
input_folder = "../input/"
output_folder = "./output/"
temp_files = "./temp/"

plotter.CoolerPlotter
pelg

def main():

    test_dataframes()

    pass




def generate_plots1(chrom_name = "chr1",
                    bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                    modle_5000_file_name = "outfile_binsize5000_tcd10.cool" , 
                    modle_1000_file_name = "outfile_binsize1000_tcd10.cool",
                    valid_5000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/5000",
                    valid_1000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/1000",
                    start = "400_000", 
                    end = "1_000_000"):

    coolerobj = cooler.Cooler(input_folder+"4DNFI9GMP2J8.mcool")
    exit()

    start = "1_000_000"
    end = "2_500_000"

    start_time = time.time()
    print("Starting function 'generate_plots1'...")

    bed_file_path = input_folder + bed_file_name
    modle_5000_file_path = input_folder + modle_5000_file_name
    modle_1000_file_path = input_folder + modle_1000_file_name
    valid_5000_file_path = input_folder + valid_5000_file_name
    valid_1000_file_path = input_folder + valid_1000_file_name


    ## * Acquire cooler-objects and their respective 2D selectors
    cooler_object_time_start = time.time()
    print("Acquiring cooler objects")
    modle_5000_cool_object = pelg.File_Functions.read_mcool_file(modle_5000_file_path)
    modle_5000_cool_2D_selector = modle_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    modle_1000_cool_object = pelg.File_Functions.read_mcool_file(modle_1000_file_path)
    modle_1000_cool_2D_selector = modle_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    valid_5000_cool_object = pelg.File_Functions.read_mcool_file(valid_5000_file_path)
    valid_5000_cool_2D_selector = valid_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    valid_1000_cool_object = pelg.File_Functions.read_mcool_file(valid_1000_file_path)
    valid_1000_cool_2D_selector = valid_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)


    cooler_object_time_end = time.time() - cooler_object_time_start
    print("Finished acquiring cooler objects in:", cooler_object_time_end,"s")

    modle_5000_numpyfilename = f'{modle_5000_file_name}.{chrom_name}.{start}.{end}'
    modle_5000_matrix_raw_numpyfilename = f'{modle_5000_numpyfilename}.matrix.raw.npy'
    modle_5000_matrix_re_numpyfilename = f'{modle_5000_numpyfilename}.matrix.re.npy'
    modle_5000_dataframe_raw_numpyfilename = f'{modle_5000_numpyfilename}.dataframe.raw.npy'
    modle_5000_dataframe_re_numpyfilename = f'{modle_5000_numpyfilename}.dataframe.re.npy'
    modle_5000_distance_dataframe_raw_numpyfilename = f'{modle_5000_numpyfilename}.distancedataframe.raw.npy'
    modle_5000_distance_dataframe_re_numpyfilename = f'{modle_5000_numpyfilename}.distancedataframe.re.npy'

    modle_1000_numpyfilename = f'{modle_1000_file_name}.{chrom_name}.{start}.{end}'
    modle_1000_matrix_raw_numpyfilename = f'{modle_1000_numpyfilename}.matrix.raw.npy'
    modle_1000_matrix_re_numpyfilename = f'{modle_1000_numpyfilename}.matrix.re.npy'
    modle_1000_dataframe_raw_numpyfilename = f'{modle_1000_numpyfilename}.dataframe.raw.npy'
    modle_1000_dataframe_re_numpyfilename = f'{modle_1000_numpyfilename}.dataframe.re.npy'
    modle_1000_scaled_5000_dataframe_re_numpyfilename = f'{modle_1000_numpyfilename}.scaled.5000.dataframe.re.npy'
    modle_1000_distance_dataframe_raw_numpyfilename = f'{modle_1000_numpyfilename}.distancedataframe.raw.npy'
    modle_1000_distance_dataframe_re_numpyfilename = f'{modle_1000_numpyfilename}.distancedataframe.re.npy'

    valid_5000_numpyfilename = f'{valid_5000_file_name}.{chrom_name}.{start}.{end}'.replace('::/resolutions/','')
    valid_5000_matrix_raw_numpyfilename = f'{valid_5000_numpyfilename}.matrix.raw.npy'
    valid_5000_matrix_re_numpyfilename = f'{valid_5000_numpyfilename}.matrix.re.npy'
    valid_5000_dataframe_raw_numpyfilename = f'{valid_5000_numpyfilename}.dataframe.raw.npy'
    valid_5000_dataframe_re_numpyfilename = f'{valid_5000_numpyfilename}.dataframe.re.npy'
    valid_5000_distance_dataframe_raw_numpyfilename = f'{valid_5000_numpyfilename}.distancedataframe.raw.npy'
    valid_5000_distance_dataframe_re_numpyfilename = f'{valid_5000_numpyfilename}.distancedataframe.re.npy'

    valid_1000_numpyfilename = f'{valid_1000_file_name}.{chrom_name}.{start}.{end}'.replace('::/resolutions/','')
    valid_1000_matrix_raw_numpyfilename = f'{valid_1000_numpyfilename}.matrix.raw.npy'
    valid_1000_matrix_re_numpyfilename = f'{valid_1000_numpyfilename}.matrix.re.npy'
    valid_1000_dataframe_raw_numpyfilename = f'{valid_1000_numpyfilename}.dataframe.raw.npy'
    valid_1000_dataframe_re_numpyfilename = f'{valid_1000_numpyfilename}.dataframe.re.npy'
    valid_1000_scaled_5000_dataframe_re_numpyfilename = f'{valid_1000_numpyfilename}.scaled.5000.dataframe.re.npy'
    valid_1000_distance_dataframe_raw_numpyfilename = f'{valid_1000_numpyfilename}.distancedataframe.raw.npy'
    valid_1000_distance_dataframe_re_numpyfilename = f'{valid_1000_numpyfilename}.distancedataframe.re.npy'

    print("Fetching dataframes...")
    dataframe_start_time = time.time()

    ## *  MoDLE 5000
    raw_modle_5000_dataframes = load_or_create_dataframe(raw_dataframe_filepath=modle_5000_dataframe_raw_numpyfilename,
                                                        re_dataframe_filepath=modle_5000_dataframe_re_numpyfilename,
                                                        raw_distance_dataframe_filepath=modle_5000_distance_dataframe_raw_numpyfilename,
                                                        re_distance_dataframe_filepath=modle_5000_distance_dataframe_re_numpyfilename,
                                                        cooler_file_path=modle_5000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end)

    raw_modle_5000_dataframe = raw_modle_5000_dataframes[0]
    re_modle_5000_dataframe = raw_modle_5000_dataframes[1]
    raw_modle_5000_distance_dataframe_expanded = raw_modle_5000_dataframes[2]
    re_modle_5000_distance_dataframe_expanded = raw_modle_5000_dataframes[3]

    ## * MoDLE 1000
    raw_modle_1000_dataframes = load_or_create_dataframe(raw_dataframe_filepath=modle_1000_dataframe_raw_numpyfilename,
                                                        re_dataframe_filepath=modle_1000_dataframe_re_numpyfilename,
                                                        raw_distance_dataframe_filepath=modle_1000_distance_dataframe_raw_numpyfilename,
                                                        re_distance_dataframe_filepath=modle_1000_distance_dataframe_re_numpyfilename,
                                                        cooler_file_path=modle_1000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end,)

    raw_modle_1000_dataframe = raw_modle_1000_dataframes[0]
    re_modle_1000_dataframe = raw_modle_1000_dataframes[1]
    raw_modle_1000_distance_dataframe_expanded = raw_modle_1000_dataframes[2]
    re_modle_1000_distance_dataframe_expanded = raw_modle_1000_dataframes[3]

    ## * MoDLE 1000 scaled
    raw_modle_1000_dataframes = load_or_create_dataframe(re_dataframe_filepath=modle_1000_scaled_5000_dataframe_re_numpyfilename,
                                                        cooler_file_path=modle_1000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end,
                                                        scale=5000)

    re_modle_1000_dataframe_scaled_5000 = raw_modle_1000_dataframes[0]


    ## * Valid 5000
    raw_valid_5000_dataframes = load_or_create_dataframe(raw_dataframe_filepath=valid_5000_dataframe_raw_numpyfilename,
                                                        re_dataframe_filepath=valid_5000_dataframe_re_numpyfilename,
                                                        raw_distance_dataframe_filepath=valid_5000_distance_dataframe_raw_numpyfilename,
                                                        re_distance_dataframe_filepath=valid_5000_distance_dataframe_re_numpyfilename,
                                                        cooler_file_path=valid_5000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end)

    raw_valid_5000_dataframe = raw_valid_5000_dataframes[0]
    re_valid_5000_dataframe = raw_valid_5000_dataframes[1]
    raw_valid_5000_distance_dataframe_expanded = raw_valid_5000_dataframes[2]
    re_valid_5000_distance_dataframe_expanded = raw_valid_5000_dataframes[3]

    ## * Valid 1000
    raw_valid_1000_dataframes = load_or_create_dataframe(raw_dataframe_filepath=valid_1000_dataframe_raw_numpyfilename,
                                                        re_dataframe_filepath=valid_1000_dataframe_re_numpyfilename,
                                                        raw_distance_dataframe_filepath=valid_1000_distance_dataframe_raw_numpyfilename,
                                                        re_distance_dataframe_filepath=valid_1000_distance_dataframe_re_numpyfilename,
                                                        cooler_file_path=valid_1000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end)

    raw_valid_1000_dataframe = raw_valid_1000_dataframes[0]
    re_valid_1000_dataframe = raw_valid_1000_dataframes[1]
    raw_valid_1000_distance_dataframe_expanded = raw_valid_1000_dataframes[2]
    re_valid_1000_distance_dataframe_expanded = raw_valid_1000_dataframes[3]

    ## * MoDLE 1000 scaled
    raw_modle_1000_dataframes = load_or_create_dataframe(re_dataframe_filepath=valid_1000_scaled_5000_dataframe_re_numpyfilename,
                                                        cooler_file_path=valid_1000_file_path,
                                                        bed_file_path=bed_file_path,
                                                        chrom_name=chrom_name,
                                                        start=start,
                                                        end=end,
                                                        scale=5000)

    re_valid_1000_dataframe_scaled_5000 = raw_modle_1000_dataframes[0]

    dataframe_end_time = time.time() - dataframe_start_time
    print("Finished fetching dataframes in:",dataframe_end_time)
    
    # * Plots
    print("Plotting....")
    print("Creating distance plots")
    distance_plots_start_time = time.time()
    # Create CoolerPlotter object
    coolerplotter = plotter.CoolerPlotter()
    
    ## * Plot average distance between bins
    coolerplotter.average_distance_plot(raw_valid_5000_distance_dataframe_expanded,line_lable="Validation 5000",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(raw_valid_1000_distance_dataframe_expanded,line_lable="Validation 1000",
                                        show = False)
    coolerplotter.average_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000",
                                        show = False)
    coolerplotter.average_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000",
                                        plotname="Average counts based distance between bins")

    
    # ## * Plot average distance only between bins containing regulatory elements
    coolerplotter.average_distance_plot(re_valid_5000_distance_dataframe_expanded,line_lable="Validation 5000",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(re_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000",
                                        show = True, plotname="Average counts based distance between 5000 resolution bins with regulatory elements")
    coolerplotter.average_distance_plot(re_valid_1000_distance_dataframe_expanded,line_lable="Validation 1000",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(re_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000",
                                        show = True, plotname="Average counts based distance between 1000 resolution bins with regulatory elements")
                                        
    ## * Plot total counts of all distances
    coolerplotter.totalcounts_distance_plot(raw_valid_5000_distance_dataframe_expanded,line_lable="Validation 5000",  
                                            newplot=True,show = False)
    coolerplotter.totalcounts_distance_plot(raw_valid_1000_distance_dataframe_expanded,line_lable="Validation 1000",  
                                            show = False)
    coolerplotter.totalcounts_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000", 
                                            show = False)
    coolerplotter.totalcounts_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000", 
                                            plotname="Total counts based on distance between bins")
    

    ## * Plot standarddeviation of all distances
    coolerplotter.standard_deviation_distance_plot(raw_valid_5000_distance_dataframe_expanded,line_lable="Validation 5000",  
                                            newplot=True,show = False)
    coolerplotter.standard_deviation_distance_plot(raw_valid_1000_distance_dataframe_expanded,line_lable="Validation 1000",  
                                            show = False)
    coolerplotter.standard_deviation_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000", 
                                            show = False)
    coolerplotter.standard_deviation_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000", 
                                            plotname="Logarithmic standarddeviation based on distance between bins.")

    ## * Plot standarddeviation of all distances
    coolerplotter.standard_deviation_distance_plot(raw_valid_5000_distance_dataframe_expanded,line_lable="Validation 5000",  
                                            newplot=True,show = False,logY=False)
    coolerplotter.standard_deviation_distance_plot(raw_valid_1000_distance_dataframe_expanded,line_lable="Validation 1000",  
                                            show = False,logY=False)
    coolerplotter.standard_deviation_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000", 
                                            show = False,logY=False)
    coolerplotter.standard_deviation_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000", 
                                            plotname="Standarddeviation based on distance between bins.",logY=False)


    distance_plots_end_time = time.time() - distance_plots_start_time
    print("Finished distance plots in:", distance_plots_end_time)


    ## * Plot matrices
    print("Plotting contact matrices")
    contact_matrix_plots_start_time = time.time()
    ## * Plot raw matrices
    if start and end:
        modle_5000_annotation="MoDLE:5000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        modle_1000_annotation="MoDLE:1000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        valid_5000_annotation="H1:5000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        valid_1000_annotation="H1:1000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
    else:
        modle_5000_annotation="MoDLE:5000 " + chrom_name
        modle_1000_annotation="MoDLE:1000 " + chrom_name
        valid_5000_annotation="H1:5000 " + chrom_name
        valid_1000_annotation="H1:1000 " + chrom_name

    coolerplotter.simple_contact_matrix_plot(in_data = modle_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_5000_annotation,
                                        out_filepath="matrix_MoDLE_5000_tcd10_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = modle_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_1000_annotation,
                                        out_filepath="matrix_MoDLE_1000_tcd5_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = valid_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_5000_annotation,
                                        out_filepath="matrix_H1_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = valid_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_1000_annotation,
                                        out_filepath="matrix_H1_1000_chr19_400KB-1MB.png",open_in_viewer=True)
    
    ## * Plot regulatory matrices


    # re_modle_1000_dataframe_scaled_5000 = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(re_modle_1000_dataframe,5000)
    # re_valid_1000_dataframe_scaled_5000 = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(re_valid_1000_dataframe,5000)
    modle_1000_annotation_scaled = modle_1000_annotation + ' scaled to 5000'
    valid_1000_annotation_scaled = valid_1000_annotation + ' scaled to 5000'


    coolerplotter.simple_contact_matrix_plot(in_data = re_modle_5000_dataframe,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=modle_5000_annotation,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_MoDLE_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_valid_5000_dataframe,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=valid_5000_annotation,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_H1_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_modle_1000_dataframe_scaled_5000,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=modle_1000_annotation_scaled,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_MoDLE_1000_scaled_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_valid_1000_dataframe_scaled_5000,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=valid_1000_annotation_scaled,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_H1_1000_scaled_5000_chr19_400KB-1MB.png",open_in_viewer=True)


    contact_matrix_plots_end_time = time.time() - contact_matrix_plots_start_time
    print("Finished plotting contact matrices in:",contact_matrix_plots_end_time, "s")

    ## * Converting to matrix 

    re_modle_5000_matrix = coolerplotter.dataframe_to_matrix(re_modle_5000_dataframe,start,end)
    re_valid_5000_matrix = coolerplotter.dataframe_to_matrix(re_valid_5000_dataframe,start,end)
    re_modle_1000_matrix_scaled_5000 = coolerplotter.dataframe_to_matrix(re_modle_1000_dataframe_scaled_5000,start,end)
    re_valid_1000_matrix_scaled_5000 = coolerplotter.dataframe_to_matrix(re_valid_1000_dataframe_scaled_5000,start,end)


    matrices = [re_modle_5000_matrix, re_modle_1000_matrix_scaled_5000, re_valid_5000_matrix,  re_valid_1000_matrix_scaled_5000]
    
    ## * Scatterplots

    for i in range(len(matrices)):
        for j in range(i+1, len(matrices)):
            xlable = ""
            ylable = ""
            title  = ""
            out_filepath="scatterplot_"
            if i == 0: 
                out_filepath = out_filepath + "MoDLE_5000_"
                xlable="MoDLE 5000bp"
                title="MoDLE 5000bp | "
            if i == 1: 
                out_filepath = out_filepath + "MoDLE_1000_"
                xlable="MoDLE 1000bp"
                title="MoDLE 1000bp | "
            if i == 2: 
                out_filepath = out_filepath + "Valid_5000_"
                xlable="Valid 5000bp"
                title="Valid 5000bp | "
            if i == 3: 
                out_filepath = out_filepath + "Valid_1000_"
                xlable="Valid 1000bp"
                title="Valid 1000bp | "

            if j == 0: 
                out_filepath = out_filepath + "MoDLE_5000.png"
                ylable="MoDLE 5000bp"
                title=title+" MoDLE 5000bp"
            if j == 1: 
                out_filepath = out_filepath + "MoDLE_1000.png"
                ylable="MoDLE 1000bp"
                title=title+" MoDLE 1000bp"
            if j == 2: 
                out_filepath = out_filepath + "Valid_5000.png"
                ylable="Valid 5000bp"
                title=title+" Valid 5000bp"
            if j == 3: 
                out_filepath = out_filepath + "Valid_1000.png"
                ylable="Valid 1000bp"
                title=title+" Valid 1000bp"

            coolerplotter.scatter_plot(matrices[i],matrices[j],xlable,ylable,title,out_filepath,open_in_viewer=True,)


    end_time = time.time() - start_time
    print("Function 'generate_plots1' finished runtime in:", end_time, "s")

    return



def compare_tcd_plots(chrom_name = "chr19",
                    bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                    tcd1_filename = "outfile_binsize1000.cool" , 
                    tcd5_filename = "outfile_binsize1000_tcd5.cool",
                    tcd10_filename = "outfile_binsize1000_tcd10.cool",
                    valid_filename = "4DNFI9GMP2J8.mcool::/resolutions/1000",
                    start = "400_000", 
                    end = "1_000_000"):



    bed_file_path = input_folder + bed_file_name

    tcd1_filepath = input_folder + tcd1_filename
    tcd1_numpyfilename = f'{tcd1_filename}.{chrom_name}.{start}.{end}'
    tcd1_matrix_raw_numpyfilename = f'{tcd1_numpyfilename}.matrix.raw.npy'
    tcd1_matrix_re_numpyfilename = f'{tcd1_numpyfilename}.matrix.re.npy'
    tcd1_dataframe_raw_numpyfilename = f'{tcd1_numpyfilename}.dataframe.raw.npy'
    tcd1_dataframe_re_numpyfilename = f'{tcd1_numpyfilename}.dataframe.re.npy'
    tcd1_distance_dataframe_raw_numpyfilename = f'{tcd1_numpyfilename}.distancedataframe.raw.npy'
    tcd1_distance_dataframe_re_numpyfilename = f'{tcd1_numpyfilename}.distancedataframe.re.npy'

    tcd5_filepath = input_folder + tcd5_filename
    tcd5_numpyfilename = f'{tcd5_filename}.{chrom_name}.{start}.{end}'
    tcd5_matrix_raw_numpyfilename = f'{tcd5_numpyfilename}.matrix.raw.npy'
    tcd5_matrix_re_numpyfilename = f'{tcd5_numpyfilename}.matrix.re.npy'
    tcd5_dataframe_raw_numpyfilename = f'{tcd5_numpyfilename}.dataframe.raw.npy'
    tcd5_dataframe_re_numpyfilename = f'{tcd5_numpyfilename}.dataframe.re.npy'
    tcd5_distance_dataframe_raw_numpyfilename = f'{tcd5_numpyfilename}.distancedataframe.raw.npy'
    tcd5_distance_dataframe_re_numpyfilename = f'{tcd5_numpyfilename}.distancedataframe.re.npy'

    tcd10_filepath = input_folder + tcd10_filename
    tcd10_numpyfilename = f'{tcd10_filename}.{chrom_name}.{start}.{end}'
    tcd10_matrix_raw_numpyfilename = f'{tcd10_numpyfilename}.matrix.raw.npy'
    tcd10_matrix_re_numpyfilename = f'{tcd10_numpyfilename}.matrix.re.npy'
    tcd10_dataframe_raw_numpyfilename = f'{tcd10_numpyfilename}.dataframe.raw.npy'
    tcd10_dataframe_re_numpyfilename = f'{tcd10_numpyfilename}.dataframe.re.npy'
    tcd10_distance_dataframe_raw_numpyfilename = f'{tcd10_numpyfilename}.distancedataframe.raw.npy'
    tcd10_distance_dataframe_re_numpyfilename = f'{tcd10_numpyfilename}.distancedataframe.re.npy'

    valid_filepath = input_folder + valid_filename
    valid_numpyfilename = f'{valid_filename}.{chrom_name}.{start}.{end}'.replace("::/resolutions/","")
    valid_matrix_raw_numpyfilename = f'{valid_numpyfilename}.matrix.raw.npy'
    valid_matrix_re_numpyfilename = f'{valid_numpyfilename}.matrix.re.npy'
    valid_dataframe_raw_numpyfilename = f'{valid_numpyfilename}.dataframe.raw.npy'
    valid_dataframe_re_numpyfilename = f'{valid_numpyfilename}.dataframe.re.npy'
    valid_distance_dataframe_raw_numpyfilename = f'{valid_numpyfilename}.distancedataframe.raw.npy'
    valid_distance_dataframe_re_numpyfilename = f'{valid_numpyfilename}.distancedataframe.re.npy'

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


    tcd1_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd1_distance_dataframe_raw_numpyfilename,
                                        re_distance_dataframe_filepath= tcd1_distance_dataframe_re_numpyfilename,
                                        cooler_file_path=tcd1_filepath,
                                        bed_file_path=bed_file_path,
                                        chrom_name=chrom_name,
                                        start = start,
                                        end = end,)
    


    tcd5_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd5_distance_dataframe_raw_numpyfilename,
                                        re_distance_dataframe_filepath= tcd5_distance_dataframe_re_numpyfilename,
                                        cooler_file_path=tcd5_filepath,
                                        bed_file_path=bed_file_path,
                                        chrom_name=chrom_name,
                                        start = start,
                                        end = end,)
    
    tcd10_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= tcd10_distance_dataframe_raw_numpyfilename,
                                        re_distance_dataframe_filepath= tcd10_distance_dataframe_re_numpyfilename,
                                        cooler_file_path=tcd10_filepath,
                                        bed_file_path=bed_file_path,
                                        chrom_name=chrom_name,
                                        start = start,
                                        end = end,)

    valid_dataframes = load_or_create_dataframe(raw_distance_dataframe_filepath= valid_distance_dataframe_raw_numpyfilename,
                                    re_distance_dataframe_filepath= valid_distance_dataframe_re_numpyfilename,
                                    cooler_file_path=valid_filepath,
                                    bed_file_path=bed_file_path,
                                    chrom_name=chrom_name,
                                    start = start,
                                    end = end,)

    tcd1_distance_raw_dataframe = tcd1_dataframes[0]
    tcd1_distance_re_dataframe = tcd1_dataframes[1]

    tcd5_distance_raw_dataframe = tcd5_dataframes[0]
    tcd5_distance_re_dataframe = tcd5_dataframes[1]

    tcd10_distance_raw_dataframe = tcd10_dataframes[0]
    tcd10_distance_re_dataframe = tcd10_dataframes[1]

    valid_distance_raw_dataframe = valid_dataframes[0]
    valid_distance_re_dataframe = valid_dataframes[1]

    ## * Plot average distance between bins

    coolerplotter.average_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000:tcd10",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000:tcd5",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000:tcd1",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="Valid 1000",
                                        show = True, plotname="Average counts based distance between bins")
    
    ## * Plot average distance only between bins containing regulatory elements
    coolerplotter.average_distance_plot(dataframe = tcd10_distance_re_dataframe,line_lable="MoDLE 1000:tcd10",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(dataframe = tcd5_distance_re_dataframe,line_lable="MoDLE 1000:tcd5",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = tcd1_distance_re_dataframe,line_lable="MoDLE 1000:tcd1",
                                        show = False)
    coolerplotter.average_distance_plot(dataframe = valid_distance_re_dataframe,line_lable="Valid 1000",
                                        show = True, plotname="Average counts based distance between bins with regulatory elements")


    ## * Plot total counts of all distances
    coolerplotter.totalcounts_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000:tcd10",
                                        newplot=True,show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000:tcd5",
                                        show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000:tcd1",
                                        show = False)
    coolerplotter.totalcounts_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="Valid 1000",
                                        show = True, plotname="Total counts based distance between bins with regulatory elements")


    ## * Plot total counts of all distances
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000:tcd10",
                                        newplot=True,show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000:tcd5",
                                        show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000:tcd1",
                                        show = False)
    coolerplotter.standard_deviation_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="Valid 1000",
                                        show = True, plotname="Logarithmic standard deviation for regulatory elements with distance")



    coolerplotter.standard_deviation_distance_plot(dataframe = tcd10_distance_raw_dataframe,line_lable="MoDLE 1000:tcd10",
                                        newplot=True,show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd5_distance_raw_dataframe,line_lable="MoDLE 1000:tcd5",
                                        show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = tcd1_distance_raw_dataframe,line_lable="MoDLE 1000:tcd1",
                                        show = False, logY=False)
    coolerplotter.standard_deviation_distance_plot(dataframe = valid_distance_raw_dataframe,line_lable="Valid 1000",
                                        show = True, plotname="Standard deviation for regulatory elements with distance", logY=False)    

    

    tcd1_matrices = load_or_create_matrix(raw_matrix_filename=tcd1_matrix_raw_numpyfilename,
                                            re_matrix_filename=tcd1_matrix_re_numpyfilename,
                                            cooler_file_path=tcd1_filepath,
                                            chrom_name=chrom_name,
                                            start=start,
                                            end=end,
                                            bed_file_path=bed_file_path,
                                            raw_distance_dataframe_filepath=tcd1_distance_raw_dataframe,
                                            re_distance_dataframe_filepath=tcd1_distance_re_dataframe
                                        )
    
    tcd1_raw_matrix = tcd1_matrices[0]
    tcd1_re_matrix = tcd1_matrices[1]

    tcd5_matrices = load_or_create_matrix(raw_matrix_filename=tcd5_matrix_raw_numpyfilename,
                                            re_matrix_filename=tcd5_matrix_re_numpyfilename,
                                            cooler_file_path=tcd5_filepath,
                                            chrom_name=chrom_name,
                                            start=start,
                                            end=end,
                                            bed_file_path=bed_file_path,
                                            raw_distance_dataframe_filepath=tcd5_distance_raw_dataframe,
                                            re_distance_dataframe_filepath=tcd5_distance_re_dataframe
                                        )
    
    tcd5_raw_matrix = tcd5_matrices[0]
    tcd5_re_matrix = tcd5_matrices[1]

    tcd10_matrices = load_or_create_matrix(raw_matrix_filename=tcd10_matrix_raw_numpyfilename,
                                            re_matrix_filename=tcd10_matrix_re_numpyfilename,
                                            cooler_file_path=tcd10_filepath,
                                            chrom_name=chrom_name,
                                            start=start,
                                            end=end,
                                            bed_file_path=bed_file_path,
                                            raw_distance_dataframe_filepath=tcd10_distance_raw_dataframe,
                                            re_distance_dataframe_filepath=tcd10_distance_re_dataframe
                                        )
    
    tcd10_raw_matrix = tcd10_matrices[0]
    tcd10_re_matrix = tcd10_matrices[1]

    valid_matrices = load_or_create_matrix(raw_matrix_filename=valid_matrix_raw_numpyfilename,
                                            re_matrix_filename=valid_matrix_re_numpyfilename,
                                            cooler_file_path=valid_filepath,
                                            chrom_name=chrom_name,
                                            start=start,
                                            end=end,
                                            bed_file_path=bed_file_path,
                                            raw_distance_dataframe_filepath=valid_distance_raw_dataframe,
                                            re_distance_dataframe_filepath=valid_distance_re_dataframe
                                        )
    
    valid_raw_matrix = valid_matrices[0]
    valid_re_matrix = valid_matrices[1]

    tcd1_standarddeviation = tcd1_distance_raw_dataframe["standarddeviation"].to_numpy()
    tcd5_standarddeviation = tcd5_distance_raw_dataframe["standarddeviation"].to_numpy()
    tcd10_standarddeviation = tcd10_distance_raw_dataframe["standarddeviation"].to_numpy()

#    coolerplotter.scatter_plot(tcd1_standarddeviation,tcd5_standarddeviation,"tcd1","tcd5","Sandarddeviation correlation","standarddeviation_scatterplot.png",open_in_viewer=True)

    matrices = [tcd1_raw_matrix,tcd5_raw_matrix,tcd10_raw_matrix]
    for row1 in range(len(matrices)):
        for row2 in range(row1+1, len(matrices)):
            corr_coeff = np.corrcoef(matrices[row1].flatten(), matrices[row2].flatten())[0,1]
            print(f"Correlation coefficient between raw matrix {row1+1} and raw matrix {row2+1}: {corr_coeff}")

    for row in range(len(matrices)):
        corr_coeff_valid = np.corrcoef(matrices[row].flatten(), valid_raw_matrix.flatten())[0,1]
        print(f"Correlation coefficient between raw matrix {row+1} and raw valid matrix: {corr_coeff_valid}")


    matrices = [tcd1_re_matrix,tcd5_re_matrix,tcd10_re_matrix]
    for row1 in range(len(matrices)):
        for row2 in range(row1+1, len(matrices)):
            corr_coeff = np.corrcoef(matrices[row1].flatten(), matrices[row2].flatten())[0,1]
            print(f"Correlation coefficient between re matrix {row1+1} and re matrix {row2+1}: {corr_coeff}")


    for row in range(len(matrices)):
        corr_coeff_valid = np.corrcoef(matrices[row].flatten(), valid_re_matrix.flatten())[0,1]
        print(f"Correlation coefficient between re matrix {row+1} and re valid matrix: {corr_coeff_valid}")




def load_or_create_matrix(raw_matrix_filename : str = False,
                    re_matrix_filename : str = False, 
                    raw_distance_matrix_filename : str = False,
                    re_distance_matrix_filename : str = False,
                    cooler_file_path : str = False, 
                    bed_file_path : str  = False, 
                    chrom_name  = False, 
                    start = False, 
                    end = False, 
                    scale : int = 0,
                    cache_dataframe : bool = False,
                    raw_dataframe_filepath: str = False, 
                    re_dataframe_filepath: str = False, 
                    raw_distance_dataframe_filepath: str = False,
                    re_distance_dataframe_filepath: str = False, ):

    def missing_param(param : any, param_name : str):
        if param: return False
        else: 
            print("load_matrix_set missing param:", param_name)
            return True

    return_matrices = []

    #TODO Move dataframe_to_matrix to File_Functions file
    coolerplotter = plotter.CoolerPlotter()

    if cooler_file_path: 
        cooler_object = cooler.Cooler(cooler_file_path)
        selector = cooler_object.matrix(balance=False, as_pixels=True, join=True)
        resolution = cooler_object.binsize

    re_dataframe = False

    ## * Raw matrix 
    if raw_matrix_filename:
        loaded_matrix = load_matrix(raw_matrix_filename)
        if not isinstance(loaded_matrix,np.ndarray):
            print("File",raw_matrix_filename,"not found. Creating new matrix.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): print(chrom_name); return False;         
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;  
            
            if raw_dataframe_filepath:
                loaded_dataframe = load_dataframe(raw_dataframe_filepath)    
                if not isinstance(loaded_dataframe,pd.DataFrame):
                    if start and end: raw_dataframe = selector.fetch((chrom_name,start,end))
                    else: raw_dataframe = selector.fetch((chrom_name))
                    if scale: raw_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(raw_dataframe,target_resolution = scale)
                    if cache_dataframe and raw_dataframe_filepath: save_dataframe(raw_dataframe,raw_dataframe_filepath)
                else: 
                    raw_dataframe = load_dataframe
            else:
                if start and end: raw_dataframe = selector.fetch((chrom_name,start,end))
                else: raw_dataframe = selector.fetch((chrom_name))
                if scale: raw_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(raw_dataframe,target_resolution = scale)
            
            raw_matrix = coolerplotter.dataframe_to_matrix(raw_dataframe)
            save_matrix(raw_matrix, raw_matrix_filename)
        else:                         
            raw_matrix = loaded_matrix
        return_matrices.append(raw_matrix)
    
    
    ## * Regulatory matrix
    if re_matrix_filename:
        loaded_matrix = load_matrix(re_matrix_filename)
        if not isinstance(loaded_matrix,np.ndarray):     
            print("File",re_matrix_filename,"not found. Creating new matrix.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
            #if missing_param(start, f'{start=}'.split('=')[0]): return False;  
            #if missing_param(end, f'{end=}'.split('=')[0]): return False;  
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            
            if re_dataframe_filepath: 
                loaded_dataframe = load_dataframe(re_dataframe_filepath)    
                if not isinstance(loaded_dataframe,pd.DataFrame):    
                    re_dataframe = pelg.Dataframe_Functions.filter_cooler_to_regelement_dataframe(cooler_file_path, bed_file_path, chrom_name, resolution, start, end)
                    if scale: re_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(re_dataframe,target_resolution = scale)
                    if cache_dataframe and re_dataframe_filepath: save_dataframe(re_dataframe,re_dataframe_filepath)
                else:
                    re_dataframe = loaded_dataframe
            else:
                re_dataframe = pelg.Dataframe_Functions.filter_cooler_to_regelement_dataframe(cooler_file_path, bed_file_path, chrom_name, resolution, start, end)
                if scale: re_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(re_dataframe,target_resolution = scale)
            
            re_matrix = coolerplotter.dataframe_to_matrix(re_dataframe)
            save_matrix(re_matrix, re_matrix_filename)
        else:                          
            re_matrix = loaded_matrix
        return_matrices.append(re_matrix)


    ## * Raw Distance matrix
    if raw_distance_matrix_filename:
        loaded_distancematrix = load_matrix(raw_distance_matrix_filename)
        if not isinstance(loaded_distancematrix,np.ndarray):
            print("File",raw_distance_matrix_filename,"not found. Creating new matrix.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False
            if not isinstance(raw_dataframe,pd.DataFrame): 
                if start and end: raw_dataframe = selector.fetch((chrom_name,start,end))
                else: raw_dataframe = selector.fetch((chrom_name))
            raw_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_dataframe)
            raw_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(raw_distance_dataframe)
            raw_distance_matrix = raw_distance_dataframe_expanded.values
            
            if cache_dataframe and raw_distance_dataframe_filepath: save_dataframe(raw_distance_dataframe,raw_distance_dataframe_filepath)
            save_matrix(raw_distance_matrix, raw_distance_matrix_filename)
        else:
            raw_distance_matrix = loaded_distancematrix
        return_matrices.append(raw_distance_matrix)


    ## * Regulatory Distance matrix
    if re_distance_matrix_filename:
        loaded_distancematrix = load_matrix(re_distance_matrix_filename)
        if not isinstance(loaded_distancematrix,np.ndarray):
            print("File",re_distance_matrix_filename,"not found. Creating new matrix.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
            #if missing_param(start, f'{start=}'.split('=')[0]): return False;  
            #if missing_param(end, f'{end=}'.split('=')[0]): return False;  
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            if not isinstance(re_dataframe,pd.DataFrame): re_dataframe = pelg.Dataframe_Functions.filter_cooler_to_regelement_dataframe(cooler_file_path, bed_file_path, chrom_name, resolution, start, end)

            re_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_dataframe)
            re_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(re_distance_dataframe)
            re_distance_matrix = re_distance_dataframe_expanded.values
            
            if cache_dataframe and re_distance_dataframe_filepath: save_dataframe(re_distance_dataframe,re_distance_dataframe_filepath)
            save_matrix(re_distance_matrix, re_distance_matrix_filename)
        else:
            re_distance_matrix = loaded_distancematrix
        return_matrices.append(re_distance_matrix)

    return return_matrices



def load_or_create_dataframe(raw_dataframe_filepath : str  = False, 
                            re_dataframe_filepath : str  = False, 
                            raw_distance_dataframe_filepath : str  = False, 
                            re_distance_dataframe_filepath : str  = False,
                            scale : int = 0, 
                            cooler_file_path : str = False, 
                            bed_file_path : str = False, 
                            chrom_name : str = False,
                            start : str = False,
                            end : str = False

                            ):
    """Attempts to load dataframe. If not available, will instead 

    Args:
        dataframe_filepath (str): _description_
        bedfile_path (_type_): _description_
        scale (int, optional): _description_. Defaults to 0.
    """
    
    def missing_param(param : any, param_name : str):
        if param: return False
        else: 
            print("load_matrix_set missing param:", param_name)
            return True

    #TODO Move dataframe_to_matrix to File_Functions file
    coolerplotter = plotter.CoolerPlotter()

    if cooler_file_path: 
        cooler_object = cooler.Cooler(cooler_file_path)
        selector = cooler_object.matrix(balance=False, as_pixels=True, join=True)
        resolution = cooler_object.binsize

    raw_dataframe = False
    re_dataframe = False

    returned_dataframes = []

    ## * Raw dataframe
    if raw_dataframe_filepath:
        raw_dataframe = load_dataframe(raw_dataframe_filepath)
        if not isinstance(raw_dataframe,pd.DataFrame): 
            print("Creating new dataframe.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
            if missing_param(start, f'{start=}'.split('=')[0]): return False;  
            if missing_param(end, f'{end=}'.split('=')[0]): return False;  
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            if start and end: raw_dataframe = selector.fetch((chrom_name,start,end))
            else: raw_dataframe = selector.fetch((chrom_name))
            if scale: raw_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(raw_dataframe,target_resolution = scale)
            save_dataframe(raw_dataframe, raw_dataframe_filepath)
        returned_dataframes.append(raw_dataframe)

    ## * Regulatory dataframe
    if re_dataframe_filepath:
        re_dataframe = load_dataframe(re_dataframe_filepath)
        if not isinstance(re_dataframe,pd.DataFrame):
            print("Creating new dataframe.")
            if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
            if missing_param(start, f'{start=}'.split('=')[0]): return False;  
            if missing_param(end, f'{end=}'.split('=')[0]): return False;  
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            re_dataframe = pelg.Dataframe_Functions.filter_cooler_to_regelement_dataframe(cooler_file_path, bed_file_path, chrom_name, resolution, start, end)
            if scale: re_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(re_dataframe,target_resolution = scale)
            save_dataframe(re_dataframe, re_dataframe_filepath)
        returned_dataframes.append(re_dataframe)

    ## * Raw distance dataframe
    if raw_distance_dataframe_filepath:
        raw_distance_dataframe = load_dataframe(raw_distance_dataframe_filepath)
        if not isinstance(raw_distance_dataframe,pd.DataFrame):
            print("Creating new dataframe.")
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            if not isinstance(raw_dataframe,pd.DataFrame):
                if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
                if missing_param(start, f'{start=}'.split('=')[0]): return False;  
                if missing_param(end, f'{end=}'.split('=')[0]): return False;  
                if start and end: raw_dataframe = selector.fetch((chrom_name,start,end))
                else: raw_dataframe = selector.fetch((chrom_name))
            raw_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_dataframe)
            raw_distance_dataframe = pelg.Dataframe_Functions.expand_distance_dataframe(raw_distance_dataframe)
            save_dataframe(raw_distance_dataframe, raw_distance_dataframe_filepath)
        if 'distance' in raw_distance_dataframe.columns: raw_distance_dataframe = raw_distance_dataframe.set_index('distance')
        returned_dataframes.append(raw_distance_dataframe)

    ## * Regulatory distance dataframe lol
    if re_distance_dataframe_filepath:
        re_distance_dataframe = load_dataframe(re_distance_dataframe_filepath)
        if not isinstance(re_distance_dataframe,pd.DataFrame):
            print("Creating new dataframe.")
            if missing_param(cooler_file_path, f'{cooler_file_path=}'.split('=')[0]): return False;   
            if missing_param(bed_file_path, f'{bed_file_path=}'.split('=')[0]): return False; 
            if not isinstance(re_dataframe,pd.DataFrame):
                if missing_param(chrom_name, f'{chrom_name=}'.split('=')[0]): return False;  
                if missing_param(start, f'{start=}'.split('=')[0]): return False;  
                if missing_param(end, f'{end=}'.split('=')[0]): return False; 
                re_dataframe = pelg.Dataframe_Functions.filter_cooler_to_regelement_dataframe(cooler_file_path, bed_file_path, chrom_name, resolution, start, end)
                if scale: re_dataframe = pelg.Dataframe_Functions.downscale_dataframe_to_resolution(re_dataframe,target_resolution = scale)
            re_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_dataframe)
            re_distance_dataframe = pelg.Dataframe_Functions.expand_distance_dataframe(re_distance_dataframe)
            save_dataframe(re_distance_dataframe, re_distance_dataframe_filepath)
        if 'distance' in re_distance_dataframe.columns: re_distance_dataframe = re_distance_dataframe.set_index('distance')
        returned_dataframes.append(re_distance_dataframe)
    
    return returned_dataframes

append_to_file_lock = threading.Lock()

def append_text_to_file(text : str, filename : str):
    with append_to_file_lock: 
        timenow = datetime.now()
        current_time = timenow.strftime("%H:%M:%S")

        timestamped_text = "[" + current_time + "]: " + text
        file_object = open(filename, "a")
        print(timestamped_text)
        file_object.write(timestamped_text + "\n")
        file_object.close()


if __name__ == "__main__":
    main()

