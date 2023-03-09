from generate_plots import *
import cooler

def main():
    #calculate_average_contact_density_memory_issue()
    calculate_average_contact_density(cooler_file_path=input_folder + "H1-hESC.mcool::/resolutions/50000")

def calculate_average_contact_density(cooler_file_path : str,
                            diagonal_width : int = 3_000_000):
    """Calculates the average contact density of the diagonal of a cooler file. Every bin within a {diagonal_width}
    distance of the exact diagonal, is counted as being on the diagonal. Bins are added together, then divided by the 
    total number of bins on the diagonal. 

    Can be used get an approximate value for use in MoDLE's target-contact-density paramter. 
    Keep in mind that MoDLE applies this average to every chromosome individually. But the actual average density of a 
    micro-c matrix may vary greatly.

    Args:
        cooler_file_path (str): The path to the cooler file. 
        diagonal_width (int, optional): The width of the diagonal. Defaults to 3_000_000.
    """

    

    cooler_obj = cooler.Cooler(cooler_file_path)
    print(cooler_file_path)
    print(cooler_obj)
    chrom_names = cooler_obj.chromnames
    resolution = cooler_obj.binsize
    diagonal_bin_width = diagonal_width // resolution

    contact_density_dict = {}

    for chrom_name in chrom_names:

        

        chrom_size = cooler_obj.chromsizes[chrom_name]

        fetch_start = 0
        # Fetch twice the size of the diagonal width. 
        # This is small enough to fit in memory, and big enough to ensure we don't lose any pixels, as long as we do 
        # some overlap with the next fetch. 
        fetch_end = diagonal_width + diagonal_width

        matrix = cooler_obj.matrix(balance=False).fetch(f'{chrom_name}:{fetch_start}-{fetch_end}')
        print("Chrom name:", chrom_name)
        print("size of chrom: ",chrom_size)
        print("resolution: ",resolution)
        
        diagonal_bin_width = diagonal_width // resolution
        length_of_one_row =  chrom_size // resolution
        current_start_bin = 0


        total_contacts = 0
        nonzero_bins = 0
        for row_index, row in enumerate(matrix):
            for col_index, col in enumerate(row):
                if col_index > row_index: break
                if abs(row_index - col_index) > diagonal_bin_width: continue

                bin_count = matrix[row_index][col_index]

                total_contacts += bin_count
                if bin_count > 0: nonzero_bins += 1

        

        while fetch_end < chrom_size - 1:
            fetch_start = fetch_start + diagonal_width
            fetch_end = fetch_end + diagonal_width
            if fetch_end > chrom_size: fetch_end = chrom_size-1
            print(fetch_start,fetch_end)

            matrix = cooler_obj.matrix(balance=False).fetch(f'{chrom_name}:{fetch_start}-{fetch_end}')

            current_start_bin = current_start_bin + diagonal_bin_width

            for row_index, row in enumerate(matrix):
                for col_index, col in enumerate(row):

                    if row_index <= len(matrix) / 2 and col_index <= len(matrix) / 2: continue
                    if col_index > row_index: continue
                    if abs(row_index - col_index) > diagonal_bin_width: continue

                    bin_count = matrix[row_index][col_index]


                    total_contacts += bin_count
                    if bin_count > 0: nonzero_bins += 1

        bins_on_diagonal = 0
        # Calculate the total number of bins within {diagonal_bin_width} of the diagonal
        for i in range(length_of_one_row + 1):
            # Add row and column less than diagonal_bin_width away from coordinate. 
            bins_on_row_or_column = diagonal_bin_width
            # Make sure we dont go beyond the end of the matrix. If we do, only count bins up until the edge
            if i + diagonal_bin_width > length_of_one_row: bins_on_row_or_column = length_of_one_row - i
            # Make sure we keep point exactly on diagonal
            bins_on_diagonal += bins_on_row_or_column + 1

        # We are interested in the average contact density on the diagonal
        average_on_diagonal = total_contacts / bins_on_diagonal

        contact_density_dict[chrom_name] = average_on_diagonal

    # Get the average of all 
    genome_wide_average_on_diagonal = 0
    for val in contact_density_dict.values():
        genome_wide_average_on_diagonal += val
    genome_wide_average_on_diagonal = genome_wide_average_on_diagonal / len(contact_density_dict.values())

    contact_density_dict['genome_wide'] = genome_wide_average_on_diagonal

    print(contact_density_dict)

    return contact_density_dict

def calculate_average_contact_density_memory_issue(chrom_name = "chr19",
                            cooler_file_name = "H1-hESC.mcool::/resolutions/5000",
                            diagonal_width = 3_000_000):
    """Deprecated. Do not use.

    Args:
        chrom_name (str, optional): _description_. Defaults to "chr19".
        cooler_file_name (str, optional): _description_. Defaults to "H1-hESC.mcool::/resolutions/5000".
        diagonal_width (_type_, optional): _description_. Defaults to 3_000_000.
    """
    # cooler_file_name = "outfile_binsize1000_tcd10.cool"
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
        bin_1_end = row[3]
        bin_2_start = row[5]
        bin_count = row[7]

        if abs(bin_1_end - bin_2_start) > diagonal_width: 
            print(bin_1_end,bin_2_start)
            continue

        total_contacts += bin_count
        nonzero_bins += 1

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

