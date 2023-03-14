import pandas as pd
import numpy as np
from Constants import *


def create_distance_dataframe(dataframe: pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    """Creates a dataframe where each row has the distance between two bins, and the contact frequency between them"""

    print("Creating distance dataframe.")

    format1 = False
    format2 = False

    def _format1():
        distance_matrix = pd.DataFrame(columns=["distance","count"])
        for row in dataframe.itertuples():
            bin1_start_index = dataframe.columns.get_loc("start1")
            bin1_start = row[bin1_start_index + 1]
            bin2_start_index = dataframe.columns.get_loc("start2")
            bin2_start = row[bin2_start_index + 1]

            distance = abs(bin2_start - bin1_start)
            count_index = dataframe.columns.get_loc("count")
            count = row[count_index + 1]

            distance_matrix.loc[len(distance_matrix)] = [distance,count]

        distance_matrix = distance_matrix.convert_dtypes()
        return distance_matrix

    def _format2():
        distance_matrix = pd.DataFrame(columns=["distance","count"])
        for row in dataframe.itertuples():
            bin1_start_index = dataframe.columns.get_loc("prom_start")
            bin1_start = row[bin1_start_index + 1]
            bin2_start_index = dataframe.columns.get_loc("enh_start")
            bin2_start = row[bin2_start_index + 1]

            distance = abs(bin2_start - bin1_start)
            count_index = dataframe.columns.get_loc("modle_count")
            count = row[count_index + 1]

            distance_matrix.loc[len(distance_matrix)] = [distance,count]

        distance_matrix = distance_matrix.convert_dtypes()
        return distance_matrix

    if ("start1" in dataframe.columns) and ("start2" in dataframe.columns) and ("count" in dataframe.columns):
        format1 = True
    elif ("bin1_start" in dataframe.columns) and ("bin2_start" in dataframe.columns) and ("modle_count" in dataframe.columns):
        format2 = True
    else:    
        print("Invalid format.")
    
    if format1: 
        return _format1()
    if format2:
        return _format2()

def expand_distance_dataframe(dataframe : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    """ Looks for 'distance" and 'count' columsn in the given dataframe. 
    If present, expands the dataframe to also include the foolowing columns:
            'distance', 'mincount', 'maxcount,'
            'averagecount', 'mediancount', 'standarddeviation',
            'totalcount', 'numberofcounts', 'allcounts'


    """

    print("Expanding distance dataframe with more statistics.")

    # Check if dataframe already has correct format. If not, reformat. 
    if not "distance" in dataframe.columns and not "count" in dataframe.columns:
        print("NO__")

    
    # Find all unique distances
    uqdistances = np.unique(dataframe['distance'].to_numpy())
    emptyarray = [None for i in uqdistances]
    emptyarrayarray = [[] for i in uqdistances]

    length = len(uqdistances)
    # Temporary dictionary to easily make a dataframe (because pandas is a pain)
    dictionary = {'distance': uqdistances, 'mincount': emptyarray, 'maxcount': emptyarray,
            'averagecount': emptyarray, 'mediancount': emptyarray, 'standarddeviation': emptyarray,
            'totalcount': emptyarray, 'numberofcounts': emptyarray, 'allcounts': emptyarrayarray}


    # Create a dataframe with distance as the index
    new_dataframe = pd.DataFrame.from_dict(dictionary)
    new_dataframe['distance'] = new_dataframe['distance'].astype('int64')
    new_dataframe = new_dataframe.set_index('distance')

    # Iterate through the old dataframe
    for row in dataframe.itertuples():
        distance = int(row[1])
        count = row[2]

        # Insert all counts with the same distance into an array at the same index 
        new_dataframe.at[distance,'allcounts'].append(count)
    # Iterate through new dataframe and update all data based on the counts the 'allcounts' array
    for row in new_dataframe.itertuples():
        allcounts = row[DATAFRAME_COLUMNS_STATISTICAL.index('allcounts')]
        mincount = min(allcounts)
        maxcounts = max(allcounts)
        averagecount = np.mean(allcounts)
        mediancount = np.median(allcounts)
        standarddeviation = np.std(allcounts)
        totalcount = np.sum(allcounts)
        numberofcounts = len(allcounts)

        # Disables warning from the following line
        np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 
        new_dataframe.loc[row[0]] = [mincount,maxcounts,averagecount,mediancount,standarddeviation,totalcount,numberofcounts,allcounts]

    return new_dataframe
