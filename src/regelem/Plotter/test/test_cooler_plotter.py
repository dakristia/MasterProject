from cooler_ploter import *

def test_simple_dataframe_plot():
    cooler_file_path : str = "output5.cool" 

    cooler_ploter = CoolerPloter()

    cooler_object, cooler_2Dselector = cooler_ploter.read_Modlefile(cooler_file_path)

    #cooler_ploter.simple_dataframe_plot(None,cooler_object, open_in_viewer = True, start = "75000", end = "248905000")
    cooler_ploter.simple_dataframe_plot(None,cooler_object, open_in_viewer = True, start = "75000", end = "500000", out_filepath="plot_withLog50_000_cmapfall.png")


def test_simple_plot():
    cooler_file_path : str = "output5.cool" 

    cooler_ploter = CoolerPloter()

    cooler_object, cooler_2Dselector = cooler_ploter.read_Modlefile(cooler_file_path)

    cooler_ploter.simple_plot(cooler_object, open_in_viewer = True)


def test_compare_plot():
    matrix1 = np.matrix([[30, 30], [30, 30]])
    matrix2 = np.matrix([[29, 31], [20, 40]])

    cooler_ploter = CoolerPloter()

    cooler_ploter.compare_plot(matrix1,matrix2)

def test_dataframe_to_matrix():
    test_dataframe = pd.DataFrame(Columns= DATAFRAME_COLUMNS_INTERNAL)
    
    #cooler_plotter.dataframe_to_matrix()

def main():
    #test_simple_plot()
    #test_compare_plot()
    test_simple_dataframe_plot()

main()