import numpy as np
import pandas as pd

class dataCollection:
    label : str
    chrom_name : str
    start : int 
    end : int
    raw_dataframe : pd.DataFrame
    regulatory_dataframe : pd.DataFrame
    distance_dataframe : pd.DataFrame
    raw_matrix : np.ndarray
    regulatory_matrix : np.ndarray
    distance_matrix : np.ndarray
    statistics : str #placeholder
    misc : np.ndarray 