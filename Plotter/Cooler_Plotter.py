import cooler 
import numpy as np #type:ignore
import pandas as pd #type:ignore
import matplotlib
matplotlib.use('Agg') # Make sure matplot doesn't try to display the plots
import matplotlib.pyplot as plt #type:ignore
import matplotlib.scale as scl  #type ignore
import sys
import subprocess
from matplotlib.colors import LogNorm
from Constants import *

OUTFOLDER = "./output"

class CoolerPlotter:
    def __init__(self):
        pass

    #! This should either be deleted or moved to PromoterEnhancerListGenerator module
    def compare_plot(self,matrix1 : np.matrix,matrix2):

        comparison_matrix = matrix1 - matrix2

        plt.matshow(comparison_matrix, cmap='YlOrRd')
        plt.savefig('comparison.png')

        return comparison_matrix

    #! This should either be deleted or moved to PromoterEnhancerListGenerator module
    def normalize_matrix(self, matrix : list):
        highest = np.amax(matrix)        
        matrix = matrix / highest

        return matrix

    #! This sucks 
    def equalize_and_normalize_matrices(self, matrices : list):
            
        first_matrix = matrices[0]
        first_matrix_sum = np.sum(first_matrix)
        matrix_list = [first_matrix]

        for matrix in matrices[1:]:
            ratio = first_matrix_sum / np.sum(matrix)
            matrix = matrix * ratio
            matrix_list.append(matrix)
        
        for matrix in matrix_list: print("sum",np.sum(matrix))

        for index, matrix in enumerate(matrix_list):
            highest = np.amax(matrix)
            lowest = np.amin(matrix)
            matrix_scaled = np.array([(x - lowest) / (highest - lowest) for x in matrix])
            matrix_list[index] = matrix_scaled

        for matrix in matrix_list: print("sum",np.sum(matrix))

        return tuple(matrix_list)


    def equalize_matrices(self, matrices : list):
        first_matrix = matrices[0]
        first_matrix_sum = np.sum(first_matrix)
        matrix_list = [first_matrix]

        for matrix in matrices[1:]:
            ratio = first_matrix_sum / np.sum(matrix)
            matrix = matrix * ratio
            matrix_list.append(matrix)
        
        for matrix in matrix_list: print(np.sum(matrix))

        return tuple(matrix_list)

    def coorelate_matrix(self,
                        dataframe : pd.DataFrame,
                        ):
        dataframe = dataframe.convert_dtypes()
        dataframe = dataframe[dataframe.apply(lambda row: False if row[1] == 0 else True, axis=1)]

    def average_distance_plot(self,dataframe,
                                line_lable = "", 
                                plotname = "line_plot", 
                                newplot = False, 
                                show = True,
                                logY = True):

        if newplot: print("Creating new average distance plot with plotname:",plotname, " and line_label:", line_lable)
        else: print("Adding line to average distance plot:", line_lable)

        if newplot: plt.figure()
        file_name = plotname + '.png'

        index = dataframe.index.to_numpy()
        averagecount = dataframe['averagecount'].to_numpy()

        # Normalizing matrix
        averagecount = averagecount / averagecount.sum()

        # plt.grid(visible=True, which='major', axis='both')
        plt.title(plotname)
        plt.xlabel('distance'); plt.ylabel('contact_frequency')

        # plt.xscale('log')
        if logY: plt.yscale('log')
        #plt.plot(index,averagecount)

        #if line_lable: 
        plt.plot(index,averagecount,label = line_lable)

        plt.legend()

        plt.savefig(file_name)
        
        if show: self.view_plot(file_name)

    def totalcounts_distance_plot(self, dataframe, 
                                    line_lable, 
                                    plotname = "line_plot", 
                                    newplot = False, 
                                    show = True,
                                    logY = True):

        if newplot: print("Creating new total distance plot with plotname:",plotname, " and line_label:", line_lable)
        else: print("Adding line to total distance plot:", line_lable)

        if newplot: plt.figure()
        file_name = plotname + '.png'

        index = dataframe.index.to_numpy()
        totalcount = dataframe['totalcount'].to_numpy()
        
        # Normalizing matrix
        totalcount = totalcount / totalcount.sum()
        

        plt.title(plotname)
        plt.xlabel('distance'); plt.ylabel('contact_frequency')
        #plt.xscale('log')
        if logY: plt.yscale('log')
        
        #if line_lable: 
        plt.plot(index,totalcount,label = line_lable)

        plt.legend()

        plt.savefig(file_name)

        if show: self.view_plot(file_name)

    def standard_deviation_distance_plot(self,dataframe,
                                            line_lable = "", 
                                            plotname = "line_plot", 
                                            newplot = False, 
                                            show = True,
                                            logY = True):

        if newplot: print("Creating new standard deviation distance plot with plotname:",plotname, " and line_label:", line_lable)
        else: print("Adding line to standard deviation distance plot:", line_lable)

        if newplot: plt.figure()
        file_name = plotname + '.png'

        index = dataframe.index.to_numpy()
        standarddeviation = dataframe['standarddeviation'].to_numpy()

        # Normalizing matrix
        standarddeviation = standarddeviation / standarddeviation.sum()

        # plt.grid(visible=True, which='major', axis='both')
        plt.title(plotname)
        plt.xlabel('distance'); plt.ylabel('standarddeviation')
        plt.yticks([])

        #* We dont really want log scale for this 
        #plt.xscale('log')
        if logY: plt.yscale('log')

        plt.plot(index,standarddeviation,label = line_lable)

        plt.legend()

        plt.savefig(file_name)
        
        if show: self.view_plot(file_name)    


    def dataframe_to_matrix(self, 
                            dataframe : pd.DataFrame,
                            start : str = None,
                            end : str = None) -> np.ndarray:
        
        column_set = False
        
        ## ! Tacky, bad solution
        if 'Unnamed: 0' in dataframe.columns:
            dataframe = dataframe.drop(columns='Unnamed: 0')

        try:
            if (dataframe.columns == DATAFRAME_COLUMNS_COOLER).all() and not column_set:
                columns = DATAFRAME_COLUMNS_COOLER
                bin1_start_label = "start1"
                bin1_end_label = "end1"
                bin2_start_label = "start2"
                bin2_end_label = "end2"
                count_label = "count"
                column_set = True
        except ValueError:
            pass
        try:
            if (dataframe.columns == DATAFRAME_COLUMNS_INTERNAL).all() or (dataframe.columns == DATAFRAME_COLUMNS_INTERNAL).all() and not column_set:
                columns = DATAFRAME_COLUMNS_INTERNAL
                bin1_start_label = "bin1_start"
                bin1_end_label = "bin1_end"
                bin2_start_label = "bin2_start"
                bin2_end_label = "bin2_end"
                count_label = "modle_count"
                column_set = True
        except ValueError:
            pass
        if not column_set: 
            print("Function not implemented for these dataframe columns")
            return False

        resolution = abs(dataframe[bin2_start_label][0] - dataframe[bin2_end_label][0])
        if start == None: start = min(dataframe[bin2_start_label])
        else: start = int(start)

        if end == None: end = max(dataframe[bin2_end_label])
        else: end = int(end)


        length = int((end-start)/resolution)

        matrix = np.zeros(shape=[length,length])

        for row in dataframe.itertuples():
            row_index = row[0]
            bin1_start_index = columns.index(bin1_start_label)
            bin1_start = row[bin1_start_index + 1]
            bin1_end_index = columns.index(bin1_end_label)

            bin2_start_index = columns.index(bin2_start_label)
            bin2_start = row[bin2_start_index+ 1]
            bin2_end_index = columns.index(bin2_end_label)
            bin2_end = row[bin2_end_index + 1]

            matrix_row : int = int((bin1_start - start)/ resolution)
            matrix_column : int = int((bin2_start - start)/ resolution)

            modle_count_index = columns.index(count_label)
            modle_count = row[modle_count_index + 1]

            matrix[matrix_row][matrix_column] += modle_count
            
            if matrix_row != matrix_column:
                matrix[matrix_column][matrix_row] += modle_count


        return matrix

    

    def scatter_plot(self, 
                    matrix1 : np.ndarray ,
                    matrix2 : np.ndarray ,
                    xlabel : str = "",
                    ylabel : str = "",
                    title : str = "",
                    out_filepath : str = 'plot.png',
                    open_in_viewer : bool = False,) -> None:

        print("Creating scatterplot with title:", title, " xlabel:", xlabel, "ylabel", ylabel)

        matrix1 = self.normalize_matrix(matrix1)
        matrix2 = self.normalize_matrix(matrix2)


        print(matrix1[0:10],matrix2[0:10])

        while len(matrix1) > len(matrix2):
            matrix2 = np.append(matrix2,0)
        while len(matrix2) > len(matrix1):
            matrix1 = np.append(matrix1,0)




        px = 1/plt.rcParams['figure.dpi']
        plt.subplots(figsize=(600*px,600*px))
        plt.grid(visible=True, which='major', axis='both')
        plt.xlabel(xlabel); plt.ylabel(ylabel)
        plt.scatter(matrix1.flatten(),matrix2.flatten())
        plt.title(title)

        ax = plt.gca()
        plt.xscale('log')
        plt.yscale('log')        

        plt.savefig(out_filepath)
        if open_in_viewer: self.view_plot(out_filepath)

    def simple_logScale(self,
                        counts : np.ndarray,
                        out_filepath : str = "logPlot.png",
                        open_in_viewer : bool = False):
    
        if len(counts.shape) == 2: counts = counts.flatten()

        #plt.subplot(222)
        plt.plot(counts)
        plt.yscale('log')
        plt.title('log')
        plt.grid(True)

        plt.show()

        plt.savefig(out_filepath)

        self.view_plot(out_filepath)
            

    def simple_matrix_plot(self, 
                    matrix : np.ndarray ,
                    out_filepath : str = 'plot.png',
                    open_in_viewer : bool = False,
                    normalize = True,
                    axis_start = None,
                    axis_end = None) -> None:
        """Creates a matrix plot of a given matrix

        Args:
            matrix (np.ndarray): matrix to plot
            out_filepath (str, optional): filename to save plot image. Defaults to 'plot.png'.
            open_in_viewer (bool, optional): wether to open the plot in the OS's default image viewer. Defaults to False.
            normalize (bool, optional): wether to normalize the data or not. Defaults to True.
        """

        if axis_start == None: axis_start = 0

        ## 
        matrix_shape = matrix.shape
        x_shape = matrix_shape[0]

        # Set the number of ticks
        x_ticks_num = 4

        axis_total = axis_end - axis_start
        axis_interval = round(axis_total / (x_ticks_num - 1))

        x_ticks_array = np.array([])
        x_ticks_points = np.array([])

        # 
        for i in range(x_ticks_num):
            tick = axis_total - (axis_interval * i) + axis_start
            x_ticks_array = np.append(x_ticks_array, tick)

            tick_point = x_shape - (round(x_shape / (x_ticks_num - 1)) * i)
            x_ticks_points = np.append(x_ticks_points, tick_point)
        
        x_ticks_array = np.flip(x_ticks_array).astype(int)
        x_ticks_points = np.flip(x_ticks_points).astype(int)
        

        log_matrix = np.log10(matrix)
        
        norm = LogNorm(vmax=np.amax(log_matrix))

        # out_filepath = "deleteme.png"

        im = plt.matshow(log_matrix, cmap='YlOrRd', norm=norm)
        #print("Log", scl.LogScale(1))

        plt.colorbar(im, fraction=0.046, pad=0.04, label = 'counts')
        plt.xticks(ticks = x_ticks_points, labels = x_ticks_array)
        plt.yticks(ticks = x_ticks_points, labels = x_ticks_array)


        plt.savefig(out_filepath)

        if open_in_viewer: self.view_plot(out_filepath)

        return

    def matrix_to_csv(self,
                        matrix : np.ndarray,
                        outpath : str = 'csv'):
        print("Writing matrix to csv")
        pd.DataFrame(matrix).to_csv(outpath)


    def simple_dataframe_plot(self, 
                    in_dataframe : pd.DataFrame, 
                    cooler_object : cooler.api.Cooler,
                    out_filepath : str = 'plot.png',
                    open_in_viewer : bool = False,
                    chrom_name : str = 'chr1',
                    start : str = '0',
                    end : str = '100,000') -> None:
        
        out_filepath = OUTFOLDER = "/" + out_filepath

        
        fetchstr = chrom_name + ':' + start + '-' + end
        numpy_matrix = cooler_object.matrix(balance=False).fetch(fetchstr)

        norm = LogNorm(vmax=50_000)
        #cmap = YlOrRd
        #cmap = 'fall'
        #ax = numpy_matrix[0, 0]


        im = plt.matshow(numpy_matrix, cmap='YlOrRd', norm=norm)
        plt.colorbar(im ,fraction=0.046, pad=0.04, label='counts (linear)')



        plt.savefig(out_filepath)
        if open_in_viewer: self.view_plot(out_filepath)
        

        return
        cooler_2Dselector = cooler_object.matrix(balance=False, as_pixels=True, join=True)
        dataframe = cooler_2Dselector.fetch(fetchstr)
        print(dataframe)
        #! YOU LEFT OF HERE
        DATAFRAME_COLUMNS_INTERNAL = ["chrom", "enh_start", "enh_end", "prom_start", "prom_end", 
                            "enh_name", "prom_name", "bin1_start", "bin1_end", 
                            "bin2_start", "bin2_end", "modle_count", "micro-C_count"]
        reduced_dataframe = in_dataframe[[""]]

    def define_outfilename(self, 
                    chrom_name : str = None,
                    start : str = '0',
                    end : str = '100,000'
                    ) -> str:
        outstring = "plot"
        if chrom_name: outstring += "_" + chrom_name
        outstring += "_start" + start
        outstring += "_end" + end
        outstring += ".png"
        return outstring
    
    #! RENAME
    def simple_contact_matrix_plot(self,
                    in_data, 
                    chrom_name : str = '',
                    start = None,
                    end = None,
                    title = "",
                    annotation = "",
                    diagonal_line=False,
                    grid=False,
                    out_filepath : str = 'plot.png',
                    open_in_viewer : bool = False):
        """Creates a interaction plot from a cooler-object or dataframe. Regions need to be specified with 'start' and 'end' parameters

        Args:
            cooler_object (cooler.api.Cooler): A cooler object
            out_filepath (str, optional): File to save the plot. Defaults to 'plot.png'.
            open_in_viewer (bool, optional): Wether to automatically open the outfile or not. Defaults to False.
            chrom_name (str, optional): Specify chromosome. Defaults to 'chr1'.
            start (str, optional): Start of region. Defaults to '0'.
            end (str, optional): End of region. Defaults to '100,000'.
        """

            

        f,ax = plt.subplots(figsize=(7,6))
        norm = LogNorm(vmin=1, vmax=1000)

        

        ## * If a cooler object was sent in, fetch matrix. If a dataframe was sent in, convert to matrix
        if(type(in_data) == cooler.api.Cooler):
            if start != None and end != None: matrix = in_data.matrix(balance=False).fetch((chrom_name,start,end))
            else: 
                matrix = in_data.matrix(balance=False).fetch(chrom_name)
                start = 0; end = in_data.chromsizes["chr19"]
        elif(type(in_data) == pd.DataFrame):
            if start == None or end == None:
                raise Exception("When passing a pandas Dataframe," f'{start=}'.split('=')[0], f'{end=}'.split('=')[0], "variables are required.")
            matrix = self.dataframe_to_matrix(in_data,start,end)
        else:
            raise TypeError("Wrong type for variable", f'{in_data}', " Requires either cooler.api.Cooler object or pandas.DataFrame object")


        # Create correct ticks. 
        x_ticks_array = np.array([])
        x_ticks_points = np.array([])
        # Set the number of ticks in plot
        x_ticks_num = 5; axis_total = int(end)-int(start); axis_interval = axis_total/(x_ticks_num - 1); axis_start = int(start)
        x_shape = matrix.shape[0]

        for i in range(x_ticks_num):
            tick = axis_total - (axis_interval * i) + axis_start
            x_ticks_array = np.append(x_ticks_array, tick)

            tick_point = x_shape - (round(x_shape / (x_ticks_num - 1)) * i)
            x_ticks_points = np.append(x_ticks_points, tick_point)

        x_ticks_array = np.flip(x_ticks_array).astype(int)
        x_ticks_points = np.flip(x_ticks_points).astype(int)

        # Done creating ticks


        im = ax.matshow(matrix, norm=norm, cmap='copper_r')#vmax=2500)
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        x_pos, y_pos = xlim[1], ylim[0]

        plt.colorbar(im, fraction=0.046, pad=0.04, label='raw counts')
        
        plt.xticks(ticks = x_ticks_points, labels = x_ticks_array)
        plt.yticks(ticks = x_ticks_points, labels = x_ticks_array)
        plt.title(title,y=1.08)
        ax.annotate(annotation, xy=(1, 0), xycoords='axes fraction',
            xytext=(0, -30), textcoords='offset points',
            ha='right', va='bottom',
            bbox=dict(boxstyle='square,pad=0.5', fc='grey', alpha=0.5)),
            #arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        #plt.text(x_pos,y_pos,'Some text',ha='right', va='bottom',fontsize=12)
        if grid: plt.grid(visible=True, which='major', axis='both')
        #if diagonal_line: ax.plot([1, 0], [1, 0], transform=ax.transAxes)
        # Unclear what this does
        ax.xaxis.set_label_position('top')


        plt.savefig(out_filepath)

        if open_in_viewer: self.view_plot(out_filepath)



    def old_simple_cooler_plot(self, 
                    cooler_object : cooler.api.Cooler, 
                    out_filepath : str = 'plot.png',
                    open_in_viewer : bool = False,
                    chrom_name : str = 'chr1',
                    start : str = '0',
                    end : str = '100,000'
                    ) -> None:
        """Creates a interaction plot from a cooler object. Regions need to be specified with 'start' and 'end' parameters

        Args:
            cooler_object (cooler.api.Cooler): A cooler object
            out_filepath (str, optional): File to save the plot. Defaults to 'plot.png'.
            open_in_viewer (bool, optional): Wether to automatically open the outfile or not. Defaults to False.
            chrom_name (str, optional): Specify chromosome. Defaults to 'chr1'.
            start (str, optional): Start of region. Defaults to '0'.
            end (str, optional): End of region. Defaults to '100,000'.
        """

        resolution = cooler_object.info['bin-size']
        fetchstr = chrom_name + ':' + start + '-' + end

        dataframe = cooler_object.matrix(balance=False).fetch(fetchstr)
        #mat = cooler_object.matrix(balance=False).fetch('chr1:10,000-100,000')

        # Show matrix

        
        # Plot
        ax = dataframe[0, 0]
        im = plt.matshow(dataframe, cmap='YlOrRd') # Not sure what the log does
        #? im = plt.matshow(mat, cmap='YlOrRd')
        plt.colorbar(im,fraction=0.046, pad=0.04, label='counts')
        #TODO plt.legend()
        
        plt.savefig(out_filepath)

        if open_in_viewer: self.view_plot(out_filepath)

        #? Potentially a way plot specific pixels
        # from matplotlib import cm
        # pyplot.imshow(matrix, interpolation='nearest', cmap=cm.Blues)
        # pyplot.scatter([6,8], [10,7], color='red', s=40)
        # pyplot.show()


        
        

    def view_plot(self, filepath : str):
        """Attempts to opens an image file (.jpg, .png) of a plot in the system's default image viewer.

        Args:
            filepath (str): Filepath to the image file
        """
        print("Attempting to open in file viewer...");
        imageViewerFromCommandLine = {'linux':'xdg-open',
                                    'win32':'explorer',
                                    'darwin':'open'}[sys.platform]
        subprocess.run([imageViewerFromCommandLine, filepath])

    #! This shouldnt be here
    def read_Modlefile(self, cooler_file_path: str, print_stats : bool = False) -> tuple:
        """Reads the info in a .cool file and prints it to the terminal. 

        Keyword arguments:
        cooler_file_path -- a filepath to a .cool file

        Returns:
        A tuple containing the cooler object of the file, and a 2D selector for the matrix.
        """
        # Attempt to import data with cooler
        modle_cool: cooler.api.Cooler = cooler.Cooler(cooler_file_path)
        resolution: np.int64 = modle_cool.binsize
        chromnames: list = modle_cool.chromnames
        chromsizes: pd.core.series.Series = modle_cool.chromsizes
        info: dict = modle_cool.info

        if print_stats:
            print("Cooler stats -----------------------")
            print("Resolution:", resolution, "\n\n")
            print("chromnames:", chromnames, "\n\n")
            print("chromssizes:", chromsizes, "\n\n")
            print("info:", info, "\n\n")

        # This is a 2D selector object, very light weight
        cooler_2Dselector = modle_cool.matrix(
            balance=False, as_pixels=True, join=True)

        return (modle_cool, cooler_2Dselector)


