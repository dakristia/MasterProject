import seaborn #type:ignore
import matplotlib.pyplot as plt #type:ignore
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import LogNorm
from landscapes.single_objective import ackley
import plot_utils

def plot_tripple_correlation(feature1 : list, feature2 : list, feature3 : list, 
                            label1 : str = "Feature1", 
                            label2 : str = "Feature2", 
                            label3 : str = "Feature3",
                            file_name : str = "tripple_plot.png",
                            show : bool = True):

    
    data = pd.DataFrame({label1: feature1,label2: feature2,label3: feature3})

    matrix = np.empty(shape=(len(feature1),len(feature2)))

    unique_feature1 = np.unique(label1)
    unique_feature2 = np.unique(label2)

    pivot_table = data.pivot_table(values=label3, index=label1, columns=label2, fill_value=0)

    # index_counter = 0
    # for i, e in enumerate(unique_feature1):
    #     for y, e2 in enumerate(unique_feature2):
    #         matrix[i][y] = feature3[index_counter]
    #         index_counter += 1


    plt.figure(figsize=(12,14))
    #,annot_kws={'size': 16}
    seaborn.heatmap(pivot_table, annot=True, cmap='autumn', 
                    linewidths=1, 
                    norm=LogNorm(vmin=pivot_table.min().min(), 
                    vmax=pivot_table.max().max()),)

    # Set style for the plot
    #seaborn.set(style="whitegrid")



    #seaborn.pairplot(data)

    plt.savefig(file_name)

    if show: plot_utils.view_plot(file_name)

    


