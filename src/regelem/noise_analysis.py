import pandas as pd
import numpy as np
import cooler
import Constants
import files
import utils
import plot_utils
import matplotlib.pyplot as plt #type:ignore



def register_noise(reg_dataframe_path : str, cooler_path : str, resolution : int, output_path : str = "./output/dataframes/noise_dataframe.csv"):

    #TODO: Either make multiprocess or reduce the amount of cooler.fetch calls that are made.

    reg_dataframe = files.load_dataframe(reg_dataframe_path)

    cooler_object = cooler.Cooler(cooler_path)
    selector = cooler_object.matrix(balance=False)

    temp_dict = {"chrom":[],"prom_name":[],"enh_name":[],
                "bin1_start":[],"bin1_end":[],
                "bin2_start":[],"bin2_end":[],
                "count":[],"noise":[],
                "left":[],"upleft":[],"up":[], "upright":[],
                "right":[],"downright":[],"down":[],"downleft":[],
                "adjacent":[]}

    for row in reg_dataframe.itertuples():
        chrom_index = reg_dataframe.columns.get_loc("chrom") + 1
        prom_name_index = reg_dataframe.columns.get_loc("prom_name") + 1
        enh_name_index = reg_dataframe.columns.get_loc("enh_name") + 1
        bin1_start_index = reg_dataframe.columns.get_loc("bin1_start") + 1
        bin1_end_index = reg_dataframe.columns.get_loc("bin1_end") + 1
        bin2_start_index = reg_dataframe.columns.get_loc("bin2_start") + 1
        bin2_end_index = reg_dataframe.columns.get_loc("bin2_end") + 1
        modle_count_index = reg_dataframe.columns.get_loc("modle_count") + 1

        chrom = row[chrom_index]
        prom_name = row[prom_name_index]
        enh_name = row[enh_name_index]
        bin1_start = row[bin1_start_index]
        bin1_end = row[bin1_end_index]
        bin2_start = row[bin2_start_index]
        bin2_end = row[bin2_end_index]
        count = row[modle_count_index]

        fetch1_start = bin1_start - resolution
        fetch1_end = bin1_end + resolution
        fetch2_start = bin2_start - resolution
        fetch2_end = bin2_end + resolution

        string1 = f"{chrom}:{fetch1_start}-{fetch1_end}"
        string2 = f"{chrom}:{fetch2_start}-{fetch2_end}"

        # * This line is accountable for 99% of used time
        submatrix = selector.fetch(string1,string2)

        noise = np.sum(submatrix) - count

        left = upleft = up = upright = right = downright = down = downleft = None

        try: left = submatrix[1,0] 
        except IndexError: pass
        try: upleft = submatrix[0,0] 
        except IndexError: pass
        try: up = submatrix[0,1] 
        except IndexError: pass
        try: upright = submatrix[0,2] 
        except IndexError: pass
        try: right = submatrix[1,2] 
        except IndexError: pass
        try: downright = submatrix[2,2] 
        except IndexError: pass
        try: down = submatrix[2,1] 
        except IndexError: pass
        try: downleft = submatrix[2,0] 
        except IndexError: pass

        # * If pixel is right next to diagonal (one pixel in cardinal directions), 
        # * then it's mirrored duplicate would be counted as noise. Set this noise to 0
        distance = abs(bin1_start - bin2_start)
        if distance == resolution:
            if bin1_start > bin2_start:
                # * upright is duplicate, set to 0
                upright = None
            else:
                # * downleft is duplicate, set to 0
                downleft = None

        


        adjacent = 0
        for i in [left, upleft, up, upright, right, downright, down, downleft]:
            if i is not None:
                adjacent = adjacent + 1

        # temp_dict = {"prom_name":[],"enh_name":[],"bin1_start":[],"bin1_end":[],"bin2_start":[],"bin2_end":[],"count":[],"noise":[]}
        temp_dict["chrom"].append(chrom)
        temp_dict["prom_name"].append(prom_name)
        temp_dict["enh_name"].append(enh_name)
        temp_dict["bin1_start"].append(bin1_start)
        temp_dict["bin1_end"].append(bin1_end)
        temp_dict["bin2_start"].append(bin2_start)
        temp_dict["bin2_end"].append(bin2_end)
        temp_dict["count"].append(count)
        temp_dict["noise"].append(noise)
        temp_dict["left"].append(left)
        temp_dict["upleft"].append(upleft)
        temp_dict["up"].append(up)
        temp_dict["upright"].append(upright)
        temp_dict["right"].append(right)
        temp_dict["downright"].append(downright)
        temp_dict["down"].append(down)
        temp_dict["downleft"].append(downleft)
        temp_dict["adjacent"].append(adjacent)


    noise_dataframe = pd.DataFrame(temp_dict)

    files.save_dataframe(noise_dataframe,output_path)

    return noise_dataframe

def generate_bias(cooler_file_path, output_path):

    # Load the cooler file
    c = cooler.Cooler(cooler_file_path)

    # Perform ICE normalization
    bias_vector, converged = cooler.ice.iterative_correction(c,cis_only=True, ignore_diags=2)

    # Save the bias values to a file
    fragments = c.bins()[:]
    fragments['bias'] = bias_vector

    print(bias_vector)
    print(fragments)

    fragments.to_csv(output_path, sep='\t', columns=['chrom', 'start', 'end', 'bias'], index=False)
