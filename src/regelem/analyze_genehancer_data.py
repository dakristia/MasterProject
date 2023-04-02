import pandas as pd
import numpy as np
from collections import defaultdict
import traceback
import os
import files
import concurrent.futures



def _analyze(promoter_genehancer_dataframe : str, enhancer_genehancer_dataframe : str, resolution : int, output_path : str = False) -> pd.DataFrame:

        max_array_len = len(promoter_genehancer_dataframe) * len(enhancer_genehancer_dataframe) + 10 # * Add buffer to avoid IndexError because of bad code

        max_promoter_genehancer_id_length = _find_max_length_of_column(promoter_genehancer_dataframe, 'genehancer_id')
        promoter_genehancer_id_dtype = f'U{max_promoter_genehancer_id_length}'
        max_enhancer_genehancer_id_length = _find_max_length_of_column(enhancer_genehancer_dataframe, 'genehancer_id')
        enhancer_genehancer_id_dtype = f'U{max_enhancer_genehancer_id_length}'

        
        # * We can find enhancers related to a gene with this
        gene_enhancer_hash = make_gene_enhancer_hash(enhancer_genehancer_dataframe)
        #print(gene_enhancer_hash.keys())


        max_chrom_name_length = max(_find_max_length_of_column(promoter_genehancer_dataframe, '#chrom'),_find_max_length_of_column(enhancer_genehancer_dataframe, '#chrom'))
        chrom_name_dtype = f'U{max_chrom_name_length}'

        # * Preallocate memory for arrays
        chromosome_array = np.empty(shape=max_array_len, dtype=chrom_name_dtype)

        promoter_genehancer_id_array = np.empty(shape=max_array_len, dtype=promoter_genehancer_id_dtype)
        enhancer_genehancer_id_array = np.empty(shape=max_array_len, dtype=enhancer_genehancer_id_dtype)

        promoter_start_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
        promoter_end_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
        enhancer_start_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
        enhancer_end_position_array = np.empty(shape=max_array_len, dtype=np.uint32)

        # * Add pixel indexes to this array with folliwing syntax:
        # * pixel_index_array[index].append(np.array([value1,value2]))
        pixel_index_array = np.empty(shape=(max_array_len,), dtype=object)
        pixel_index_array[:] = [[] for _ in range(max_array_len)]

        assosiation_score_array = np.empty(shape=(max_array_len,), dtype=object)
        assosiation_score_array[:] = [[] for _ in range(max_array_len)]

        assosiated_gene_array = np.empty(shape=(max_array_len,),dtype=object)
        assosiated_gene_array[:] = [[] for _ in range(max_array_len)]


        re_pair_index = 0

        # * The last index of the previous promoter. Indexing from here will save ussom time later
        prev_prom_final_index = 0

        chromosome_name = ""
        promoter_genehancer_id = ""

        for promoter_row in promoter_genehancer_dataframe.itertuples():
            chromosome_name_index = promoter_genehancer_dataframe.columns.get_loc('#chrom') + 1

            # * Notfiy user whenever we start at a new chromosome
            if promoter_row[chromosome_name_index] != chromosome_name:
                print(f"Now analyzing genehancer data for chrom {promoter_row[chromosome_name_index]}")
            chromosome_name = promoter_row[chromosome_name_index]

            #!  TODO: Make _analyze assume that we only work on one chromosome?
            # * Reduce dataframe to same chromosome as promoter to reduce overall workload
            reduced_enhancer_genehancer_dataframe = enhancer_genehancer_dataframe[enhancer_genehancer_dataframe['#chrom'] == chromosome_name]

            promoter_genehancer_id_index = promoter_genehancer_dataframe.columns.get_loc('genehancer_id') + 1
            if promoter_genehancer_id != promoter_row[promoter_genehancer_id_index]:
                prev_prom_final_index = re_pair_index
            promoter_genehancer_id = promoter_row[promoter_genehancer_id_index]

            # * Get all assosiated genes of promoter
            connected_gene_column_index = promoter_genehancer_dataframe.columns.get_loc('connected_gene_id') + 1
            assosiation_score_column_index = promoter_genehancer_dataframe.columns.get_loc('connected_gene_score') + 1

            connected_gene_array = promoter_row[connected_gene_column_index]
            connected_assosiation_score_array = promoter_row[assosiation_score_column_index]

            # * Fetch promoter data
            promoter_start_loc_column_index = promoter_genehancer_dataframe.columns.get_loc('start') + 1
            promoter_end_loc_column_index = promoter_genehancer_dataframe.columns.get_loc('end') + 1

            promoter_start_loc = promoter_row[promoter_start_loc_column_index]
            promoter_end_loc = promoter_row[promoter_end_loc_column_index]

            promoter_start_bin = promoter_start_loc // resolution
            promoter_end_bin = promoter_end_loc // resolution

            # * Iterate through all assosiated genes
            for connected_gene_id, promoter_assosiation_score in zip(connected_gene_array, connected_assosiation_score_array):

                
                
                # * Find all enhancers that are also assosiated with the gene
                assosiated_enhancers = gene_enhancer_hash[connected_gene_id] # reduced_enhancer_genehancer_dataframe[reduced_enhancer_genehancer_dataframe['connected_gene_id'].apply(lambda x: connected_gene_id in x)]

                if not len(assosiated_enhancers):
                        continue

                for enhancer_genehancer_id in assosiated_enhancers:
                    duplicate_mask = reduced_enhancer_genehancer_dataframe['genehancer_id'] == enhancer_genehancer_id
                    row = reduced_enhancer_genehancer_dataframe.loc[duplicate_mask].iloc[0]
                    enhancer_start_loc = row['start']
                    enhancer_end_loc = row['end']
                    enhancer_start_bin_index = enhancer_start_loc // resolution
                    enhancer_end_bin_index = enhancer_end_loc // resolution

                    enhancer_assosiation_score = row['connected_gene_score'][np.where(row['connected_gene_id'] == connected_gene_id)][0]
                    
                    assosiation_score = (promoter_assosiation_score, enhancer_assosiation_score)


                    combined_genehancer_id_array = np.core.defchararray.add(promoter_genehancer_id_array, enhancer_genehancer_id_array)
                    duplicate_mask = combined_genehancer_id_array == (promoter_genehancer_id + enhancer_genehancer_id)

                    if duplicate_mask.any(): 
                        # * Find the duplicate index
                        i = np.where(combined_genehancer_id_array == promoter_genehancer_id + enhancer_genehancer_id)[0][0]

                        assosiated_gene_array[i].append(connected_gene_id)
                        assosiation_score_array[i].append(assosiation_score)

                        continue
                    else:
                        # * We have found promoter-enhancer pair. Add everything except pixels indexes to arrays.
                        chromosome_array[re_pair_index] = chromosome_name
                        promoter_genehancer_id_array[re_pair_index] = promoter_genehancer_id
                        enhancer_genehancer_id_array[re_pair_index] = enhancer_genehancer_id
                        promoter_start_position_array[re_pair_index] = promoter_start_loc
                        promoter_end_position_array[re_pair_index] = promoter_end_loc
                        enhancer_start_position_array[re_pair_index] = enhancer_start_loc
                        enhancer_end_position_array[re_pair_index] = enhancer_end_loc

                        assosiated_gene_array[re_pair_index].append(connected_gene_id)
                        assosiation_score_array[re_pair_index].append(assosiation_score)

                        # * Note bins that the promoters and enhancers are located in
                        for promoter_bin_index in range(promoter_start_bin, promoter_end_bin + 1):
                            for enhancer_bin_index in range(enhancer_start_bin_index, enhancer_end_bin_index + 1):
                                pixel_index = (promoter_bin_index,enhancer_bin_index)
                                pixel_index_array[re_pair_index].append(np.array([pixel_index[0],pixel_index[1]]))

                                # * Here we have every promoter-enhancer pair, their names, their start and end bins, 
                        
                        # * Done registering for pair
                        re_pair_index += 1


        # * Resize array to drop empty elements
        chromosome_array.resize((re_pair_index,))
        promoter_genehancer_id_array.resize((re_pair_index,))
        enhancer_genehancer_id_array.resize((re_pair_index,))
        promoter_start_position_array.resize((re_pair_index,))
        promoter_end_position_array.resize((re_pair_index,))
        enhancer_start_position_array.resize((re_pair_index,))
        enhancer_end_position_array.resize((re_pair_index,))
        pixel_index_array.resize((re_pair_index,))
        assosiated_gene_array.resize((re_pair_index,))
        assosiation_score_array.resize((re_pair_index,))

        temp_dictionary = {'chrom':chromosome_array, 
        'promoter_genehancer_id':promoter_genehancer_id_array, 'enhancer_genehancer_id':enhancer_genehancer_id_array,
        'promoter_start':promoter_start_position_array,'promoter_end':promoter_end_position_array,
        'enhancer_start':enhancer_start_position_array,'enhancer_end':enhancer_end_position_array,
        'pixels':pixel_index_array,
        'assosiated_genes':assosiated_gene_array, 'assosiation_scores':assosiation_score_array}

        new_dataframe = pd.DataFrame(temp_dictionary)

        if output_path: files.save_dataframe(new_dataframe,output_path,numpy_columns=['pixels','assosiated_genes','assosiation_scores'])

        return new_dataframe

def genehancer_data_to_bins_multiprocess(input_path : str, resolution : str, output_path : str = False, workers : int = 5, parts_per_chrom : int = 10, target_chrom : str = False):

    if os.path.splitext(output_path)[1] != ".csv":   
        print("Adding .csv to filepath")
        output_path = output_path + ".csv"
        print("New output path:", output_path)

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        
        temp_output_path = os.path.splitext(output_path)[0] +"/temp/"

        # * Read and reformat dataframe to make it easier to work with
        genehancer_dataframe = read_genehancer_to_df(input_path)


        promoter_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Promoter'].reset_index(drop=True)
        enhancer_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Enhancer'].reset_index(drop=True)

        #! What do I do with this? When not even they are sure if it is a Promoter or Enhancer. Ignore for now.
        promoter_enhancer_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Promoter/Enhancer'].reset_index(drop=True)

        print("Percentage of Promoter/Enhancer in dataset:", len(promoter_enhancer_genehancer_dataframe) / len(genehancer_dataframe))


        chrom_names = np.unique(genehancer_dataframe['#chrom']) 
        futures = []

        final_dataframe : pd.DataFrame = False

        for chrom_name in np.sort(chrom_names):
            chrom_promoter_genehancer_dataframe = promoter_genehancer_dataframe[promoter_genehancer_dataframe['#chrom'] == chrom_name]
            chrom_enhancer_genehancer_dataframe = enhancer_genehancer_dataframe[enhancer_genehancer_dataframe['#chrom'] == chrom_name]
            chrom_promoter_enhancer_genehancer_dataframe = promoter_enhancer_genehancer_dataframe[promoter_enhancer_genehancer_dataframe['#chrom'] == chrom_name]

            #! TODO: Remove
            if target_chrom and chrom_name.lower() != target_chrom.lower():
                continue

            total_parts = parts_per_chrom

            if total_parts > len(chrom_promoter_genehancer_dataframe):
                total_parts = len(chrom_promoter_genehancer_dataframe)
            
            len_per_part = len(chrom_promoter_genehancer_dataframe) // total_parts

            for part in range(total_parts):
                start_of_part : int = part * len_per_part
                end_of_part : int = (part + 1) * len_per_part
                temp_output_path_of_part = temp_output_path + f"{chrom_name}/{part}.csv"

                # * If last part, fetch the rest of the dataframe
                if part == total_parts-1:
                    end_of_part = len(chrom_promoter_genehancer_dataframe)

                part_chrom_promoter_genehancer_dataframe = chrom_promoter_genehancer_dataframe[start_of_part:end_of_part]

                futures.append(executor.submit(_analyze, part_chrom_promoter_genehancer_dataframe, chrom_enhancer_genehancer_dataframe, resolution, temp_output_path_of_part))

            concurrent.futures.wait(futures)

            for f in futures:
                try:
                    returned_df = f.result()

                    if not isinstance(final_dataframe,pd.DataFrame):
                        final_dataframe = returned_df

                    else:
                        final_dataframe = pd.concat([final_dataframe,returned_df],ignore_index=True)

                except Exception as e:
                    print(f"Error occurred: {e}\n{traceback.format_exc()}")

            futures = []




        final_dataframe = final_dataframe.sort_values(by=['chrom','promoter_start','enhancer_start']).reset_index(drop=True)

        print("final:",final_dataframe)

        files.save_dataframe(final_dataframe,file_path=output_path,numpy_columns=['pixels','assosiated_genes','assosiation_scores'])

        return final_dataframe


def _find_max_length_of_column(dataframe : str, column_name : str, verbose : bool = False):
    
    max_length = 0

    max_length = dataframe[column_name].astype(str).map(len).max()

    if verbose: print(f"Max length of {column_name}: {max_length}")

    return max_length

def read_genehancer_to_df(file_path : str):

    df = pd.read_csv(file_path, delimiter='\t')

    # * Make a genehancer id column
    genehancer_id_column_name = 'genehancer_id'
    genehancer_id_column_array = np.empty(shape = (len(df)),dtype='U40')

    # * Make a column that contains lists of connected gene names
    connected_gene_id_column_name = 'connected_gene_id'
    connected_gene_id_column_array = np.empty(shape = (len(df)),dtype=object)

    # * Make a column that contains lists of connected gene scores
    connected_gene_score_column_name = 'connected_gene_score'
    connected_gene_score_column_array = np.empty(shape = (len(df)),dtype=object)


    for column_index, l in enumerate(df['attributes']):
        split_line = l.split(";")

        # * Add genehancer id to column
        genehancer_id = split_line[0].split('=')[1]
        genehancer_id_column_array[column_index] = genehancer_id

        split_line = split_line[1:]

        # * Make inner list as a single entry to dataframe cell
        gene_id_array = np.empty(shape=len(split_line) // 2, dtype='U40')
        gene_score_array = np.empty(shape=len(split_line) // 2, dtype=float)

        for pos_counter, connected_gene in enumerate(split_line):
            connected_gene_split = connected_gene.split('=')

            if pos_counter % 2 == 0:
                # * Connected gene id
                gene_id = connected_gene_split[1]

                # * Integer divide pos_counter by half to get current position in array
                gene_id_array[pos_counter // 2] = gene_id
            else:
                # * Connected gene score
                gene_score = connected_gene_split[1]

                # * Integer divide pos_counter by half to get current position in array
                gene_score_array[pos_counter // 2] = gene_score

        connected_gene_id_column_array[column_index] = gene_id_array
        connected_gene_score_column_array[column_index] = gene_score_array

    df[genehancer_id_column_name] = genehancer_id_column_array
    df[connected_gene_id_column_name] = connected_gene_id_column_array
    df[connected_gene_score_column_name] = connected_gene_score_column_array

    df.drop('attributes',axis=1,inplace=True)

    return df

def make_gene_enhancer_hash(dataframe : str) -> defaultdict:
    """Create a dictionary with genes as key and related enhancers as value

    Args:
        dataframe (str): _description_

    Returns:
        defaultdict: _description_
    """
    
    # * Gene ID as key, list of enhancers as value

    gene_enhancer_dict = defaultdict(lambda: np.array([]))

    enhancer_dataframe = dataframe[dataframe['feature name'] == 'Enhancer']

    for row in enhancer_dataframe.itertuples():
        enhancer_id_index = enhancer_dataframe.columns.get_loc('genehancer_id') + 1
        connected_gene_id_index = enhancer_dataframe.columns.get_loc('connected_gene_id') + 1

        enhancer_id = row[enhancer_id_index]
        connected_genes = row[connected_gene_id_index]

        for gene_id in connected_genes:
            gene_enhancer_dict[gene_id] = np.append(gene_enhancer_dict[gene_id],enhancer_id)

    return gene_enhancer_dict

def gather_temp_files(directory_path : str):

    temp_directory_path = os.path.join(directory_path, 'temp')
    
    chrom_folders = [os.path.join(temp_directory_path, f) for f in os.listdir(temp_directory_path) if os.path.isdir(os.path.join(temp_directory_path, f))]

    #print(chrom_folders)

    for chrom_f in chrom_folders:
            #print(os.listdir(chrom_f))
            #! TODO: Change the ugliest one liner in the world
            part_dataframes = [files.load_dataframe(os.path.join(chrom_f, f), numpy_columns=['pixels','assosiated_genes','assosiation_scores'], allow_pickle=True) for f in os.listdir(chrom_f) if not os.path.isdir(os.path.join(chrom_f, f))]
            print(part_dataframes)

            chrom_dataframe = pd.concat(part_dataframes,ignore_index=True)

            print(chrom_dataframe)

            chrom_dataframe = chrom_dataframe.sort_values(by=['promoter_start','promoter_genehancer_id'])


            print(len(chrom_dataframe))
            # for part_f in part_folders:
            #     print(os.listdir(part_f))
            exit()