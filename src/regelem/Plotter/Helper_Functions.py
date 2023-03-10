import numpy as np
import pandas as pd
import cooler

def matrix_to_csv(self,
                        matrix : np.ndarray,
                        outpath : str = 'csv'):
        print("Writing matrix to csv")
        pd.DataFrame(matrix).to_csv(outpath)


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