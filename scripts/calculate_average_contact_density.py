from generate_plots import *
import cooler

def main():
    #calculate_average_contact_density_memory_issue()
    calculate_average_contact_density()

def calculate_average_contact_density(cooler_file_name = "H1-hESC.mcool::/resolutions/1000",
                            diagonal_width = 3_000_000):
    cooler_file_name = "outfile_binsize5000_tcd7.31.cool"
    cooler_file_path = input_folder + cooler_file_name

    

    cooler_obj = cooler.Cooler(cooler_file_path)
    print(cooler_file_name)
    print(cooler_obj)
    chrom_names = cooler_obj.chromnames
    resolution = cooler_obj.binsize
    diagonal_bin_width = diagonal_width // resolution

    contact_density_dict = {}

    for chrom_name in chrom_names:


        chrom_size = cooler_obj.chromsizes[chrom_name]


        current_start = 0
        current_end = diagonal_width

        fetch_start = 0
        fetch_end = diagonal_width + diagonal_width

        matrix = cooler_obj.matrix(balance=False).fetch(f'{chrom_name}:{fetch_start}-{fetch_end}')
        print(matrix)
        print(matrix.shape)
        print("size of chrom",chrom_size)
        print("resolution",resolution)
        
        end_bin = chrom_size // resolution
        diagonal_bin_width = diagonal_width // resolution
        #selector = cooler_obj.matrix(balance=False, as_pixels=True, join=True)
        #dataframe = selector.fetch([chrom_name,0,diagonal_width])
        #bins = cooler_obj.bins()[:]


        #matrix = cooler_obj.matrix(balance=True)[:]
        #print(dataframe)
        length_of_one_row =  chrom_size // resolution
        total_bins_in_chrom = length_of_one_row * length_of_one_row


        current_start_bin = 0
        current_end_bin = diagonal_width // resolution

        indexes = np.empty(shape=(0,2), dtype=int)

        total_contacts = 0
        nonzero_bins = 0
        for row_index, row in enumerate(matrix):
            for col_index, col in enumerate(row):
                if col_index > row_index: break
                if abs(row_index - col_index) > diagonal_bin_width: continue

                #indexes = np.append(indexes,np.array([[row_index,col_index]]),axis=0)
                bin_count = matrix[row_index][col_index]


                total_contacts += bin_count
                if bin_count > 0: nonzero_bins += 1


        # print(f"total:{total_contacts}")
        # print(f"nonzero_bins:{nonzero_bins}")
        # print(f"average:{total_contacts / total_bins_in_chrom}")
        # print(f"average of nonzero bins:{total_contacts / nonzero_bins}")

        # print(current_end_bin)
        # print(diagonal_bin_width)
        bins_on_diagonal = 0

        while fetch_end < chrom_size - 1:
            fetch_start = fetch_start + diagonal_width
            fetch_end = fetch_end + diagonal_width
            if fetch_end > chrom_size: fetch_end = chrom_size-1
            print(fetch_start,fetch_end)

            matrix = cooler_obj.matrix(balance=False).fetch(f'{chrom_name}:{fetch_start}-{fetch_end}')

            current_start = fetch_start
            current_start_bin = current_start_bin + diagonal_bin_width
            current_end = current_end + diagonal_width
            current_end_bin = current_end // resolution
            for row_index, row in enumerate(matrix):
                for col_index, col in enumerate(row):
                    #print(col)
                    #print(row_index,col_index)
                    #print(abs(row_index - col_index))
                    if row_index <= len(matrix) / 2 and col_index <= len(matrix) / 2: continue
                    if col_index > row_index: continue
                    #print("3")
                    if abs(row_index - col_index) > diagonal_bin_width: continue
                    #print(row_index,col_index)
                    #indexes = np.append(indexes,np.array([row_index,col_index]))
                    bin_count = matrix[row_index][col_index]


                    total_contacts += bin_count
                    if bin_count > 0: nonzero_bins += 1


        
        for i in range(length_of_one_row + 1):
            # Add row and column less than diagonal_bin_width away from coordinate. 
            bins_on_row_or_column = diagonal_bin_width
            # Make sure we dont go beyond the end of the matrix
            if i + diagonal_bin_width > length_of_one_row: bins_on_row_or_column = length_of_one_row - i
            # Make sure we keep point on diagonal
            bins_on_diagonal += bins_on_row_or_column + 1

        print(f"total:{total_contacts}")
        print(f"nonzero_bins:{nonzero_bins}")
        print(f"average:{total_contacts / total_bins_in_chrom}")
        print(f"average of nonzero bins:{total_contacts / nonzero_bins}")

        print("bins_on_diagonal:",bins_on_diagonal)

        average_on_diagonal = total_contacts / bins_on_diagonal

        print(f"average of bins on diagonal:{average_on_diagonal}")

        contact_density_dict[chrom_name] = average_on_diagonal

        print(len(indexes))
        uniq = np.unique(indexes)

        print(len(uniq))

    chr_wide_average_on_diagonal = 0
    for val in contact_density_dict.values():
        chr_wide_average_on_diagonal += val
    
    chr_wide_average_on_diagonal = chr_wide_average_on_diagonal / len(contact_density_dict.values())
    contact_density_dict['chr_wide'] = chr_wide_average_on_diagonal

    print(contact_density_dict)

def calculate_average_contact_density_memory_issue(chrom_name = "chr19",
                            cooler_file_name = "H1-hESC.mcool::/resolutions/5000",
                            diagonal_width = 3_000_000):
    cooler_file_name = "outfile_binsize1000_tcd10.cool"
    cooler_file_path = input_folder + cooler_file_name


    cooler_obj = cooler.Cooler(cooler_file_path)
    resolution = cooler_obj.binsize
    chrom_size = cooler_obj.chromsizes[chrom_name]
    print("size of chrom",chrom_size)
    diagonal_bin_width = diagonal_width / resolution

    selector = cooler_obj.matrix(balance=False, as_pixels=True, join=True)
    dataframe = selector.fetch(chrom_name)

    #matrix = cooler_obj.matrix(balance=True)[:]
    print(dataframe)
    length_of_one_row =  chrom_size // resolution
    total_bins_in_chrom = length_of_one_row * length_of_one_row
    len_of_dataframe = len(dataframe)

    print("resolution",resolution)
    print("len_of_dataframe:",len_of_dataframe)
    print("total_bins_in_chrom:",total_bins_in_chrom)

    (chrom_size / resolution) * (diagonal_width / resolution)

    
    total_contacts = 0
    nonzero_bins = 0
    for row in dataframe.itertuples():
        #num = row[]
        bin_1_end = row[3]
        bin_2_start = row[5]
        bin_count = row[7]

        

        if abs(bin_1_end - bin_2_start) > diagonal_width: 
            print(bin_1_end,bin_2_start)
            continue

        total_contacts += bin_count
        nonzero_bins += 1





    print(f"total:{total_contacts}")
    print(f"nonzero_bins:{nonzero_bins}")
    print(f"average:{total_contacts / (len_of_dataframe)}")
    print(f"average:{total_contacts / total_bins_in_chrom}")
    print(f"average of nonzero bins:{total_contacts / nonzero_bins}")
        
    end_bin = chrom_size // resolution
    diagonal_bin_width = diagonal_width / resolution
    print(diagonal_bin_width)

    bins_on_diagonal = 0

    for i in range(end_bin + 1):
        # Add row and column less than diagonal_bin_width away from coordinate. 
        bins_on_row_or_column = diagonal_bin_width
        # Make sure we dont go beyond the end of the matrix
        if i + diagonal_bin_width > end_bin: bins_on_row_or_column = end_bin - i

        # Make sure we keep point on diagonal
        bins_on_diagonal += bins_on_row_or_column + 1

    print("bins_on_diagonal:",bins_on_diagonal)
    print(f"average of bins on diagonal:{total_contacts / bins_on_diagonal}")

if __name__ == "__main__":
    main()

