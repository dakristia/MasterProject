#Generate list of enhancer-promoter pairs (in GM12878 cells) with their predicted contact frequeciesn from the MoDLE output

import warnings
import cooler                   #type:ignore
import matplotlib as mpl        #type:ignore
#mpl.use('TkAgg')    
import matplotlib.pyplot as plt #type:ignore
import numpy as np              #type:ignore
import pandas as pd             #type:ignore
from collections import defaultdict
from .File_Functions import *
# from Constants import *
# from File_Functions import *
# from General_Functions import *
#from IPython.display import display
import re
import time

# TODO: Optimize runtime. Its already pretty slow as of 18.09.2022. 
# TODO: Add type checking for functions 

# ! DO a code cleanup before you're done for the day
# Helper function. Consider deleting

# !Unfinished
def find_highest_regulatory_size(dataframe : pd.core.frame.DataFrame) -> None:
    start = dataframe["enh_start"]

    end = dataframe["enh_end"]

    #print(start)
    #print(end)
    #print(max(abs(start - end)))

    exit()
