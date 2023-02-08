import sys
from os.path import dirname, basename, isfile, join

sys.path.insert(0, '..')
from Constants import *
import numpy as np
import pandas as pd
import PromoterEnhancerListGenerator as pelg
import Plotter as plotter
import matplotlib.pyplot as plt
import cooler
import time

test_folder = "./testInput/"
input_folder = "../input/"



def main():
    generate_plots1()
    #read_extrution_barrier_file(start=400_000,end=1_000_000)



def generate_plots1(chrom_name = "chr19",
                            bed_file_name = "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed", 
                            modle_5000_file_name = "outfile_binsize5000.cool", 
                            modle_1000_file_name = "outfile_binsize1000.cool",
                            validation_5000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/5000",
                            validation_1000_file_name = "4DNFI9GMP2J8.mcool::/resolutions/1000",
                            start = "400_000", 
                            end = "1_000_000"):

    start_time = time.time()
    print("Starting function 'generate_plots1'...")

    bed_file_path = input_folder + bed_file_name
    modle_5000_file_path = input_folder + modle_5000_file_name
    modle_1000_file_path = input_folder + modle_1000_file_name
    validation_5000_file_path = input_folder + validation_5000_file_name
    validation_1000_file_path = input_folder + validation_1000_file_name

    ## * Acquire cooler-objects and their respective 2D selectors
    cooler_object_time_start = time.time()
    print("Acquiring cooler objects")
    modle_5000_cool_object = pelg.File_Functions.read_mcool_file(modle_5000_file_path)
    modle_5000_cool_2D_selector = modle_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    modle_1000_cool_object = pelg.File_Functions.read_mcool_file(modle_1000_file_path)
    modle_1000_cool_2D_selector = modle_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    validation_5000_cool_object = pelg.File_Functions.read_mcool_file(validation_5000_file_path)
    validation_5000_cool_2D_selector = validation_5000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    validation_1000_cool_object = pelg.File_Functions.read_mcool_file(validation_1000_file_path)
    validation_1000_cool_2D_selector = validation_1000_cool_object.matrix(balance=False, as_pixels=True, join=True)

    cooler_object_time_end = time.time() - cooler_object_time_start
    print("Finished acquiring cooler objects in:", cooler_object_time_end,"s")

    print("Fetching raw dataframes")
    raw_dataframe_start_time = time.time()
    ## * Raw modle data res 5000
    raw_modle_5000_dataframe = modle_5000_cool_2D_selector.fetch((chrom_name,start,end))
    raw_modle_5000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_modle_5000_dataframe)
    raw_modle_5000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(raw_modle_5000_distance_dataframe)
    
    ## * Raw modle data res 1000
    raw_modle_1000_dataframe = modle_1000_cool_2D_selector.fetch((chrom_name,start,end))
    #re_modle_1000_dataframe_collected = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(raw_modle_1000_dataframe,target_resolution = 5000)
    raw_modle_1000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_modle_1000_dataframe)
    raw_modle_1000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(raw_modle_1000_distance_dataframe)

    ## * Raw validation data res 5000
    raw_validation_5000_dataframe = validation_5000_cool_2D_selector.fetch((chrom_name,start,end))  
    raw_validation_5000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_validation_5000_dataframe)
    raw_validation_5000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(raw_validation_5000_distance_dataframe)
    
    ## * Raw validation data res 1000
    raw_validation_1000_dataframe = validation_1000_cool_2D_selector.fetch((chrom_name,start,end))  
    raw_validation_1000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(raw_validation_1000_dataframe)
    raw_validation_1000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(raw_validation_1000_distance_dataframe)

    raw_dataframe_end_time = time.time() - raw_dataframe_start_time
    print("Finished fetching raw dataframes in:",raw_dataframe_end_time,"s")

  

    re_dataframe_start_time = time.time()
    print("Fetching regulatory element dataframe in:",re_dataframe_start_time)
    ## * Regulatory elements modle data res 5000
    re_modle_5000_dataframe = pelg.Dataframe_Functions.create_regElem_df_of_cool(modle_5000_file_path, bed_file_path, chrom_name, 5000, start, end)
    re_modle_5000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_modle_5000_dataframe)
    re_modle_5000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(re_modle_5000_distance_dataframe)

    ## * Regulatory elements modle data res 1000
    re_modle_1000_dataframe = pelg.Dataframe_Functions.create_regElem_df_of_cool(modle_1000_file_path, bed_file_path, chrom_name, 1000, start, end)
    re_modle_1000_dataframe_collected_to_5000 = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(re_modle_1000_dataframe,target_resolution = 5000)
    re_modle_1000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_modle_1000_dataframe)
    re_modle_1000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(re_modle_1000_distance_dataframe)

    ## * Regulatory elements validation data res 5000
    re_validation_5000_dataframe = pelg.Dataframe_Functions.create_regElem_df_of_cool(validation_5000_file_path, bed_file_path, chrom_name, 5000, start, end)
    re_validation_5000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_validation_5000_dataframe)
    re_validation_5000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(re_validation_5000_distance_dataframe)
    
    ## * Regulatory elements validation data res 1000
    re_validation_1000_dataframe = pelg.Dataframe_Functions.create_regElem_df_of_cool(validation_1000_file_path, bed_file_path, chrom_name, 1000, start, end)
    re_validation_1000_distance_dataframe = pelg.Dataframe_Functions.create_distance_dataframe(re_validation_1000_dataframe)
    re_validation_1000_distance_dataframe_expanded = pelg.Dataframe_Functions.expand_distance_dataframe(re_validation_1000_distance_dataframe)
    
    re_dataframe_end_time = time.time() - re_dataframe_start_time
    print("Finished fetching regulatory element dataframes in:",re_dataframe_end_time)
    
    # * Plots
    print("Plotting....")
    print("Creating distance plots")
    distance_plots_start_time = time.time()
    # Create CoolerPlotter object
    coolerplotter = plotter.CoolerPlotter()
    
    ## * Plot average distance between bins
    coolerplotter.average_distance_plot(raw_validation_5000_distance_dataframe_expanded,line_lable="Validation 5000",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(raw_validation_1000_distance_dataframe_expanded,line_lable="Validation 1000",
                                        show = False)
    coolerplotter.average_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000",
                                        show = False)
    coolerplotter.average_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000",
                                        plotname="Average counts based distance between bins")
    
    # ## * Plot average distance only between bins containing regulatory elements
    coolerplotter.average_distance_plot(re_validation_5000_distance_dataframe_expanded,line_lable="Validation 5000",
                                        newplot=True,show = False)
    coolerplotter.average_distance_plot(re_validation_1000_distance_dataframe_expanded,line_lable="Validation 1000",
                                        show = False)
    coolerplotter.average_distance_plot(re_modle_5000_distance_dataframe_expanded,line_lable="MoDLE 5000",
                                        show = False)
    coolerplotter.average_distance_plot(re_modle_1000_distance_dataframe_expanded,line_lable="MoDLE 1000",
                                        plotname="Average counts based distance between bins with regulatory elements")
                                        
    ## * Plot total counts of all distances
    coolerplotter.totalcounts_distance_plot(raw_validation_5000_distance_dataframe_expanded,line_lable="Validation 5000",  
                                            newplot=True,show = False)
    coolerplotter.totalcounts_distance_plot(raw_validation_1000_distance_dataframe_expanded,line_lable="Validation 1000",  
                                            show = False)
    coolerplotter.totalcounts_distance_plot(raw_modle_5000_distance_dataframe_expanded,line_lable="MoDLE, res 5000", 
                                            show = False)
    coolerplotter.totalcounts_distance_plot(raw_modle_1000_distance_dataframe_expanded,line_lable="MoDLE, res 1000", 
                                            plotname="Total counts based on distance between bins.")
    

    distance_plots_end_time = time.time() - distance_plots_start_time
    print("Finished distance plots in:", distance_plots_end_time)
    ## * Plot matrices
    print("Plotting contact matrices")
    contact_matrix_plots_start_time = time.time()
    ## * Plot raw matrices
    if start and end:
        modle_5000_annotation="MoDLE:5000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        modle_1000_annotation="MoDLE:1000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        valid_5000_annotation="hg38:5000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        valid_1000_annotation="hg38:1000 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
    else:
        modle_5000_annotation="MoDLE:5000 " + chrom_name
        modle_1000_annotation="MoDLE:1000 " + chrom_name
        valid_5000_annotation="hg38:5000 " + chrom_name
        valid_1000_annotation="hg38:1000 " + chrom_name

    coolerplotter.simple_contact_matrix_plot(in_data = modle_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_5000_annotation,
                                        out_filepath="matrix_MoDLE_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = modle_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_1000_annotation,
                                        out_filepath="matrix_MoDLE_1000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = validation_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_5000_annotation,
                                        out_filepath="matrix_hg38_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = validation_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_1000_annotation,
                                        out_filepath="matrix_hg38_1000_chr19_400KB-1MB.png",open_in_viewer=True)

    ## * With grids enabled
    coolerplotter.simple_contact_matrix_plot(in_data = modle_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_5000_annotation,grid=True,
                                        out_filepath="matrix_MoDLE_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = modle_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=modle_1000_annotation,grid=True,
                                        out_filepath="matrix_MoDLE_1000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = validation_5000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_5000_annotation,grid=True,
                                        out_filepath="matrix_hg38_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = validation_1000_cool_object,chrom_name=chrom_name,start=start,end=end,
                                        title="Raw contacts",annotation=valid_1000_annotation,grid=True,
                                        out_filepath="matrix_hg38_1000_chr19_400KB-1MB.png",open_in_viewer=True)

    
    ## * Plot regulatory matrices


    re_modle_1000_dataframe_scaled_5000 = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(re_modle_1000_dataframe,5000)
    re_validation_1000_dataframe_scaled_5000 = pelg.Dataframe_Functions.collect_dataframe_to_match_resolution(re_validation_1000_dataframe,5000)
    modle_1000_annotation_scaled = modle_1000_annotation + ' scaled to 5000'
    valid_1000_annotation_scaled = valid_1000_annotation + ' scaled to 5000'

    print(re_modle_1000_dataframe_scaled_5000)

    coolerplotter.simple_contact_matrix_plot(in_data = re_modle_5000_dataframe,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=modle_5000_annotation,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_MoDLE_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    # coolerplotter.simple_contact_matrix_plot(in_data = re_modle_1000_dataframe,start=start,end=end,
    #                                     title="Promoter/enhancer contacts",annotation=modle_1000_annotation,diagonal_line=True,grid=True,
    #                                     out_filepath="REmatrix_MoDLE_1000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_validation_5000_dataframe,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=valid_5000_annotation,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_hg38_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    # coolerplotter.simple_contact_matrix_plot(in_data = re_validation_1000_dataframe,start=start,end=end,
    #                                     title="Promoter/enhancer contacts",annotation=valid_1000_annotation,diagonal_line=True,grid=True,
    #                                     out_filepath="REmatrix_hg38_1000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_modle_1000_dataframe_scaled_5000,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=modle_1000_annotation_scaled,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_MoDLE_1000_scaled_5000_chr19_400KB-1MB.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = re_validation_1000_dataframe_scaled_5000,start=start,end=end,
                                        title="Promoter/enhancer contacts",annotation=valid_1000_annotation_scaled,diagonal_line=True,grid=True,
                                        out_filepath="REmatrix_hg38_1000_scaled_5000_chr19_400KB-1MB.png",open_in_viewer=True)


    contact_matrix_plots_end_time = time.time() - contact_matrix_plots_start_time
    print("Finished plotting contact matrices in:",contact_matrix_plots_end_time, "s")


    ## * Correlation coeffecient

    matrix1 = coolerplotter.dataframe_to_matrix(re_modle_5000_dataframe,start,end)
    matrix2 = coolerplotter.dataframe_to_matrix(re_validation_5000_dataframe,start,end)
    matrix3 = coolerplotter.dataframe_to_matrix(re_modle_1000_dataframe_scaled_5000,start,end)
    matrix4 = coolerplotter.dataframe_to_matrix(re_validation_1000_dataframe_scaled_5000,start,end)
    # matrix5 = np.array(...)
    # matrix6 = np.array(...)

    matrices = [matrix1, matrix2, matrix3, matrix4]

    for i in range(len(matrices)):
        for j in range(i+1, len(matrices)):
            corr_coeff = np.corrcoef(matrices[i].flatten(), matrices[j].flatten())[0,1]
            print(f"Correlation coefficient between matrix {i+1} and matrix {j+1}: {corr_coeff}")


    f,ax = plt.subplots(figsize=(7,6))
    plt.title("Correlation Matrix")


    matrices = [x.flatten() for x in matrices]
    corr_matrix = np.corrcoef(matrices, rowvar=False)
    print(corr_matrix)
    plt.imshow(corr_matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()

    plt.savefig("corr_image.png")
    plotter.view_plot("corr_image.png")


    end_time = time.time() - start_time
    print("Function 'generate_plots1' finished runtime in:", end_time, "s")

    return


def generate_plotFigure_for_TADs_section_of_essay():
    #bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"
    hic_data_file_path = input_folder + "outfile_binsize5000.cool"
    valid_data_file_path = input_folder + "4DNFI9GMP2J8.mcool::/resolutions/5000"
    #hic_data_file_path = input_folder + "OUTPUT.cool"

    chrom_name = "chr19"
    start = "400_000"
    end = "1_000_000"
    #start = None
    #end = None
    #start = "0"
    #end = "1_000_000"
    #dataframe = pelg.Dataframe_Functions.extract_original_dataframe_from_cool(hic_data_file_path, chrom_name, start, end)

    cooler_object : cooler.api.Cooler = cooler.Cooler(hic_data_file_path)
    valid_cooler_object : cooler.api.Cooler = cooler.Cooler(valid_data_file_path)
    
    #send dataframe to plot function.
    coolerplotter = plotter.CoolerPlotter()

    if start and end:
        modle_annotation="MoDLE " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
        valid_annotation="hg38 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
    else:
        modle_annotation="MoDLE " + chrom_name
        valid_annotation="hg38 " + chrom_name

    coolerplotter.simple_contact_matrix_plot(in_data = cooler_object,chrom_name = chrom_name,start=start,end=end,title="MoDLE counts",annotation=modle_annotation,out_filepath="figure_TAD_section_essay.png",open_in_viewer=True)
    coolerplotter.simple_contact_matrix_plot(in_data = valid_cooler_object,chrom_name = chrom_name,start=start,end=end,title="validation counts",annotation=valid_annotation,out_filepath="figure_valid_TAD_section_essay.png",open_in_viewer=True)
#plot_distance_line_plot("chr19")

# * FINISHED
def read_extrution_barrier_file(chrom_name="chr19", 
                                resolution=1000,
                                start=None,
                                end=None,
                                extrution_barrier_bed_file="chr19_extrution_barriers.bed",
                                chrom_sizes_file = "hg38.chrom.sizes"):
    extrution_barrier_bed_file_path = input_folder + extrution_barrier_bed_file
    chrom_size_file_path = input_folder + chrom_sizes_file

    extrution_barrier_df : pd.core.frame.DataFrame = pd.read_csv(extrution_barrier_bed_file_path, 
                                                delim_whitespace=True, 
                                                index_col=False,
                                                header = None, 
                                                names = ["chrom","chromStart","chromEnd","name","score","strand"])



    chrom_size_dataframe : pd.core.frame.DataFrame = pd.read_csv(chrom_size_file_path,
                                                    delim_whitespace=True,
                                                    index_col="chrom",
                                                    header = None,
                                                    names = ["chrom", "size"])


    extrution_barrier_dict = extrution_barrier_df.to_dict()#.pop('chrom').pop('name').pop('score')
    del extrution_barrier_dict['chrom']; del extrution_barrier_dict['name']; del extrution_barrier_dict['score']

   
    chrom_size = chrom_size_dataframe.at['chr19','size']


    rows = chrom_size // resolution + 1

    array = np.zeros((rows,))
    array = np.full(rows,"",dtype=str)

    

    for i in range(0,len(extrution_barrier_dict['strand'])):
        chromStart = extrution_barrier_dict['chromStart'][i]
        chromEnd = extrution_barrier_dict['chromEnd'][i]
        strand = extrution_barrier_dict['strand'][i]

        chrom_start_bin = chromStart // resolution
        chrom_end_bin = chromEnd // resolution

        if strand not in array[chrom_start_bin]: array[chrom_start_bin] += strand
        if strand not in array[chrom_end_bin]: array[chrom_end_bin] += strand
        
    array = np.char.replace(array, '+-','b')
    array = np.char.replace(array,'-+','b')

    start_index = start // resolution
    end_index = end // resolution

    array = array[start_index:end_index+1]

    import matplotlib.pyplot as plt
    import sys
    import subprocess

    img = np.zeros((1, len(array), 3), dtype=np.uint8)

    # Set the color of each pixel based on the value of the corresponding array element
    img[0,np.where(array == '')] = [128, 128, 128]
    img[0,np.where(array == '+')] = [255, 0, 0]
    img[0,np.where(array == '-')] = [0, 0, 255]

    # Plot the image
    fig, ax = plt.subplots()
    ax.imshow(img, aspect=10, extent=(0, len(array), 0, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


    plt.savefig('./testlineplot.png')

    imageViewerFromCommandLine = {'linux':'xdg-open',
                                    'win32':'explorer',
                                    'darwin':'open'}[sys.platform]
    subprocess.run([imageViewerFromCommandLine, './testlineplot.png'])

def corr_coef():

    matrix1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    matrix2 = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])

    correlation_coefficient = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0,1]
    print("Correlation coefficient:", correlation_coefficient)

def generate_TAD_figure_for_TADs_section_of_essay(chrom_name="chr19", 
                                                    tad_boundary_file="chr19_TAD_boundary.bed"):
    tad_boundary_file_path = input_folder + tad_boundary_file

    df : pd.core.frame.DataFrame = pd.read_csv(tad_boundary_file_path, 
                                                delim_whitespace = True, 
                                                index_col = False,
                                                header = None, 
                                                names = ["chrom","chromStart","chromEnd","strength","score"])

    
    
    print(df)




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

    # print("\n\n\n\nDataframe RAW: ")
    # print(modle_tcd1_dataframe)
    # print(modle_tcd10_dataframe)
    # print(mcool_1000_dataframe)
    # print(mcool_5000_dataframe)


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



if __name__ == "__main__":
    main()
