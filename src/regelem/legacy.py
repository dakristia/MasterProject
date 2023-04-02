# def genehancer_data_to_bins(input_path : str, resolution : str, output_path : str = False):

#     genehancer_dataframe = read_genehancer_to_df(input_path)

#     enhancer_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Enhancer']
#     promoter_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Promoter']
#     promoter_enhancer_genehancer_dataframe = genehancer_dataframe[genehancer_dataframe['feature name'] == 'Promoter/Enhancer']

#     max_array_len = len(promoter_genehancer_dataframe * enhancer_genehancer_dataframe)

#     max_feature_name_length = _find_max_length_of_column(genehancer_dataframe, 'feature name')
#     feature_name_dtype = f'U{max_feature_name_length}'

#     max_chrom_name_length = _find_max_length_of_column(genehancer_dataframe, '#chrom')
#     chrom_name_dtype = f'U{max_chrom_name_length}'

#     # * Preallocate memory for arrays
#     chromosome_array = np.empty(shape=max_array_len, dtype=chrom_name_dtype)

#     promoter_genehancer_id_array = np.empty(shape=max_array_len, dtype=feature_name_dtype)
#     enhancer_genehancer_id_array = np.empty(shape=max_array_len, dtype=feature_name_dtype)

#     promoter_start_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
#     promoter_end_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
#     enhancer_start_position_array = np.empty(shape=max_array_len, dtype=np.uint32)
#     enhancer_end_position_array = np.empty(shape=max_array_len, dtype=np.uint32)

#     # * Add pixel indexes to this array with folliwing syntax:
#     # * pixel_index_array[index].append(np.array([value1,value2]))
#     pixel_index_array = np.empty(shape=(max_array_len,), dtype=object)
#     pixel_index_array[:] = [[] for _ in range(max_array_len)]

#     assosiation_score_array = np.empty(shape=(max_array_len,), dtype=object)
#     assosiation_score_array[:] = [[] for _ in range(max_array_len)]

#     assosiated_gene_array = np.empty(shape=(max_array_len,),dtype=object)
#     assosiated_gene_array[:] = [[] for _ in range(max_array_len)]

#     re_pair_index = 0

#     # * The last index of the final promoter. Indexing from here will save ussom time later
#     prev_prom_final_index = 0

#     chromosome_name = ""
#     promoter_genehancer_id = ""

#     for promoter_row in promoter_genehancer_dataframe.itertuples():
#         chromosome_name_index = promoter_genehancer_dataframe.columns.get_loc('#chrom') + 1

#         # * Notfiy user whenever we start at a new chromosome
#         if promoter_row[chromosome_name_index] != chromosome_name:
#             print(f"Now analyzing genehancer data for chrom {promoter_row[chromosome_name_index]}")
#         chromosome_name = promoter_row[chromosome_name_index]

#         promoter_genehancer_id_index = promoter_genehancer_dataframe.columns.get_loc('genehancer_id') + 1
#         if promoter_genehancer_id != promoter_row[promoter_genehancer_id_index]:
#             prev_prom_final_index = re_pair_index
#         promoter_genehancer_id = promoter_row[promoter_genehancer_id_index]

    

#         # * Get all assosiated genes of promoter
#         connected_gene_column_index = promoter_genehancer_dataframe.columns.get_loc('connected_gene_id') + 1
#         assosiation_score_column_index = promoter_genehancer_dataframe.columns.get_loc('connected_gene_score') + 1

#         connected_gene_array = promoter_row[connected_gene_column_index]
#         connected_assosiation_score_array = promoter_row[assosiation_score_column_index]

#         # * Fetch promoter data
#         promoter_start_loc_column_index = promoter_genehancer_dataframe.columns.get_loc('start') + 1
#         promoter_end_loc_column_index = promoter_genehancer_dataframe.columns.get_loc('end') + 1

#         promoter_start_loc = promoter_row[promoter_start_loc_column_index]
#         promoter_end_loc = promoter_row[promoter_end_loc_column_index]

#         promoter_start_bin = promoter_start_loc // resolution
#         promoter_end_bin = promoter_end_loc // resolution

#         # * Iterate through all assosiated genes
#         for connected_gene_id, promoter_assosiation_score in zip(connected_gene_array, connected_assosiation_score_array):

#             # * Reduce dataframe to same chromosome as promoter to reduce overall workloads
#             reduced_enhancer_genehancer_dataframe = enhancer_genehancer_dataframe[enhancer_genehancer_dataframe['#chrom'] == chromosome_name]
            
#             # * Find all enhancers that are also assosiated with the gene
#             assosiated_enhancers = reduced_enhancer_genehancer_dataframe[reduced_enhancer_genehancer_dataframe['connected_gene_id'].apply(lambda x: connected_gene_id in x)]

#             for enhancer_row in assosiated_enhancers.itertuples():

#                 # * Fetch data on assosiated enhancer
#                 enhancer_genehancer_id_index = reduced_enhancer_genehancer_dataframe.columns.get_loc('genehancer_id') + 1
#                 enhancer_genehancer_id = enhancer_row[enhancer_genehancer_id_index]

#                 enhancer_connected_gene_column_index = assosiated_enhancers.columns.get_loc('connected_gene_id')
#                 enhancer_assosiation_score_column_index = assosiated_enhancers.columns.get_loc('connected_gene_score')

#                 enhancer_connected_gene_array = enhancer_row[enhancer_connected_gene_column_index + 1]
#                 enhancer_assosiation_score_array = enhancer_row[enhancer_assosiation_score_column_index + 1]

#                 enhancer_gene_index = np.where(enhancer_connected_gene_array == connected_gene_id)[0][0]
#                 enhancer_assosiation_score = enhancer_assosiation_score_array[enhancer_gene_index]

#                 enhancer_start_loc_column_index = reduced_enhancer_genehancer_dataframe.columns.get_loc('start') + 1
#                 enhancer_end_loc_column_index = reduced_enhancer_genehancer_dataframe.columns.get_loc('end') + 1

#                 enhancer_start_loc = enhancer_row[enhancer_start_loc_column_index]
#                 enhancer_end_loc = enhancer_row[enhancer_end_loc_column_index]

#                 enhancer_start_bin_index = enhancer_start_loc // resolution
#                 enhancer_end_bin_index = enhancer_end_loc // resolution



#                 # * Assosiation score between promoter and enhancer calculated as the average between the two
#                 assosiation_score = (promoter_assosiation_score + enhancer_assosiation_score) / 2
#                 # * Keep individual assosiation scores as well
#                 assosiation_score = (assosiation_score, promoter_assosiation_score, enhancer_assosiation_score)

                
#                 # * If promoter-enhancer pair already exists in our arrays, we instead add the assosiated gene data and move on to next pair
#                 # * Do NOT index re_pair_index
#                 duplicate_pair = False
#                 if promoter_genehancer_id in np.unique(promoter_genehancer_id_array) and enhancer_genehancer_id in np.unique(enhancer_genehancer_id_array):
#                     for i, (p, e) in enumerate(zip(promoter_genehancer_id_array[prev_prom_final_index:re_pair_index],enhancer_genehancer_id_array[prev_prom_final_index:re_pair_index])):
#                         if p == promoter_genehancer_id and e == enhancer_genehancer_id:
#                             assosiated_gene_array[i].append(connected_gene_id)
#                             assosiation_score_array[i].append(assosiation_score)
#                             duplicate_pair = True
#                             break

#                 if duplicate_pair: continue

#                 # * We have found promoter-enhancer pair. Add everything except pixels indexes to arrays.
#                 chromosome_array[re_pair_index] = chromosome_name
#                 promoter_genehancer_id_array[re_pair_index] = promoter_genehancer_id
#                 enhancer_genehancer_id_array[re_pair_index] = enhancer_genehancer_id
#                 promoter_start_position_array[re_pair_index] = promoter_start_loc
#                 promoter_end_position_array[re_pair_index] = promoter_end_loc
#                 enhancer_start_position_array[re_pair_index] = enhancer_start_loc
#                 enhancer_end_position_array[re_pair_index] = enhancer_end_loc

#                 assosiated_gene_array[re_pair_index].append(connected_gene_id)
#                 assosiation_score_array[re_pair_index].append(assosiation_score)

#                 # * Note bins that the promoters and enhancers are located in
#                 for promoter_bin_index in range(promoter_start_bin, promoter_end_bin + 1):
#                     for enhancer_bin_index in range(enhancer_start_bin_index, enhancer_end_bin_index + 1):
#                         pixel_index = (promoter_bin_index,enhancer_bin_index)
#                         pixel_index_array[re_pair_index].append(np.array([pixel_index[0],pixel_index[1]]))

#                         # * Here we have every promoter-enhancer pair, their names, their start and end bins, 
                
#                 # * Done registering for pair
#                 re_pair_index += 1

#     # * Resize array to drop empty elements
#     chromosome_array.resize((re_pair_index,))
#     promoter_genehancer_id_array.resize((re_pair_index,))
#     enhancer_genehancer_id_array.resize((re_pair_index,))
#     promoter_start_position_array.resize((re_pair_index,))
#     promoter_end_position_array.resize((re_pair_index,))
#     enhancer_start_position_array.resize((re_pair_index,))
#     enhancer_end_position_array.resize((re_pair_index,))
#     pixel_index_array.resize((re_pair_index,))
#     assosiated_gene_array.resize((re_pair_index,))
#     assosiation_score_array.resize((re_pair_index,))

#     temp_dictionary = {'chrom':chromosome_array, 
#     'promoter_genehancer_id':promoter_genehancer_id_array, 'enhancer_genehancer_id':enhancer_genehancer_id_array,
#     'promoter_start':promoter_start_position_array,'promoter_end':promoter_end_position_array,
#     'enhancer_start':enhancer_start_position_array,'enhancer_end':enhancer_end_position_array,
#     'pixels':pixel_index_array,
#     'assosiated_genes':assosiated_gene_array, 'assosiation_scores':assosiation_score_array}

#     new_dataframe = pd.DataFrame(temp_dictionary)

#     if output_path: files.save_dataframe(new_dataframe,output_path,numpy_columns=['pixels','assosiated_genes','assosiation_scores'])

#     return new_dataframe


# # def note_interactions(  cool_file, 
# #                         promoter_enhancer_dataframe : pd.core.frame.DataFrame, 
# #                         chrom_name : str,
# #                         start : int = None,
# #                         end : int = None, 
# #                         resolution : int = False,
# #                         max_distance : int = False) -> pd.core.frame.DataFrame:
# #     """Looks through pixels in a cooler dataframe of a given resolution and finds those who include interactions
# #     between a promoter-like element and a enhancer-like element. Returns a dataframe. 


# #     Args:
# #         cooler_2Dselector (cooler.core.RangeSelector2D): a cooler matrix selector
# #         promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe with promoter/enhancer regions
# #         chrom_name (str): chromosome we want to look at
# #         start (int, optional): start of frame we want to look at. Defaults to None.
# #         end (int, optional): end of frame we want to look at. Defaults to None.
# #         resolution (int, optional): Deprecated.
# #         #TODO Make resolution dynamic like in note_validation
# #     Returns:
# #         pd.core.frame.DataFrame: Dataframe with counts of proposed promoter/enhancer interactions. Columns = DATAFRAME_COLUMNS_INTERNAL
# #     """

# #     print("\n\nNoting interactions for chrom:", chrom_name, "at resolution", resolution, "...")

# #     new_dataframe : pd.core.frame.DataFrame = pd.DataFrame(columns = DATAFRAME_COLUMNS_INTERNAL)

# #     cooler_2Dselector = cool_file.matrix(balance=False, as_pixels=True, join=True)
# #     chrom_size = cool_file.chromsizes[chrom_name]
# #     resolution = cool_file.binsize

# #     # TODO fetch numpy array and index manually. Roberto says this is faster than fetching dataframe
# #     if end:
# #         print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
# #         cooler_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
# #     else:
# #         print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
# #         print("Note: This may eat a lot of memory.")
# #         cooler_dataframe = cooler_2Dselector.fetch((chrom_name))

# #     ## *  Copy part of dataframe with relevant chrom_name 
# #     reduced_dataframe = promoter_enhancer_dataframe.loc[promoter_enhancer_dataframe["chrom"] == chrom_name]
# #     ## * split dataframe into pls and els dataframes
# #     promoter_dataframe, enhancer_dataframe = split_df_to_pls_els(reduced_dataframe)


# #     ## * Iterate through promoters and find bins of correlating enhancers
# #     print("Finding correlating promoter and enhancer bins.")
# #     for row in promoter_dataframe.itertuples():
# #         prom_start = row[2]
# #         prom_end = row[3]
# #         prom_name = row[4]
# #         score = row[5] #!Unused
# #         strand = row[6] #!Unused

# #         # Get the lower end of bins    
# #         prom_rounded_down = round_up_and_down(prom_start,resolution)[0]
# #         # Get the higher end of bins
# #         prom_rounded_up = round_up_and_down(prom_end,resolution)[1]

# #         ## * Getting all bins where the PLS is within the range of the first region of bin
# #         promoter_interaction_dataframe1 = cooler_dataframe.loc[cooler_dataframe["start1"] >= prom_rounded_down].loc[cooler_dataframe["end1"] <= prom_rounded_up]

# #         ## * Getting all bins where the PLS is within the range of the second region of bin
# #         promoter_interaction_dataframe2 = cooler_dataframe.loc[cooler_dataframe["start2"] >= prom_rounded_down].loc[cooler_dataframe["end2"] <= prom_rounded_up]


# #         ## * Iterating through all rows to find relevant enhancers that are in second region, while promoter is in first region
# #         for row2 in promoter_interaction_dataframe1.itertuples():
# #             enh_rounded_down = row2[5]
# #             enh_rounded_up = row2[6]
# #             count = row2[7]
            

# #             ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
# #             diff = 0
# #             if prom_rounded_up < enh_rounded_down: 
# #                 diff = abs(prom_rounded_up - enh_rounded_down)
                
# #             elif enh_rounded_up < prom_rounded_down:
# #                 diff = abs(enh_rounded_up - prom_rounded_down)

# #             if diff > DISTANCE_LIMIT: 
# #                 #print(f"Diff {diff} is higher than {DISTANCE_LIMIT}. Skipping row.")
# #                 continue

# #             ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
# #             #enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] >= enh_rounded_down].loc[enhancer_dataframe["chromStart"] < enh_rounded_up]
# #             enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]

# #             ## * Iterate through all relevant enhancers and insert into new dataframe
# #             for row3 in enhancer_hits_dataframe.itertuples():
# #                 enh_start = row3[2]
# #                 enh_end = row3[3]
# #                 enh_name = row3[4]
                
# #                 ## * Last count is -1 as we currently dont have data to insert
# #                 internal_columns = DATAFRAME_COLUMNS_INTERNAL

# #                 ## * Make a list out of data to insert into new dataframe
# #                 ## * has to be a list of list
# #                 #TODO: Change column order so promoters are first, and enhancers second.
# #                 input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
# #                             enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
# #                             prom_rounded_down, prom_rounded_up, count, -1]]

# #                 input_df = pd.DataFrame(input_list, columns = internal_columns)
# #                 new_dataframe = pd.concat([new_dataframe,input_df])
        
# #         ## * Iterating through all rows to find relevant enhancers that are in first region, while promoter is in second region        
# #         for row2 in promoter_interaction_dataframe2.itertuples():
# #             enh_rounded_down = row2[2]
# #             enh_rounded_up = row2[3]
# #             count = row2[7]



# #             ## * If distance between promoter/enhancer is longer than DISTANCE_LIMIT, skip
# #             diff = abs(prom_rounded_down - enh_rounded_down)
# #             if diff > DISTANCE_LIMIT: continue

# #             ## * Gettings all enhancers that are within bins that interacted with the PLS's bin
# #             enhancer_hits_dataframe = enhancer_dataframe.loc[enhancer_dataframe["chromStart"] <= enh_rounded_up].loc[enhancer_dataframe["chromEnd"] > enh_rounded_down]

# #             ## * Iterate through all relevant enhancers and insert into new dataframe
# #             for row3 in enhancer_hits_dataframe.itertuples():
# #                 enh_start = row3[2]
# #                 enh_end = row3[3]
# #                 enh_name = row3[4]

# #                 internal_columns = DATAFRAME_COLUMNS_INTERNAL

# #                 ## * Make a list out of data to insert into new dataframe
# #                 ## * has to be a list of list
# #                 #TODO: Change column order so promoters are first, and enhancers second.
# #                 input_list = [[chrom_name, enh_start, enh_end, prom_start, prom_end, 
# #                             enh_name, prom_name, enh_rounded_down, enh_rounded_up, 
# #                             prom_rounded_down, prom_rounded_up, count, -1]]

# #                 input_df = pd.DataFrame(input_list, columns = internal_columns)
# #                 new_dataframe = pd.concat([new_dataframe,input_df])


# #     print("Sorting dataframe and dropping duplicates.")
# #     new_dataframe = new_dataframe.drop_duplicates().sort_values(by=['prom_start','enh_start']).reset_index(drop=True)

# #     print(new_dataframe)
# #     return new_dataframe

# # def note_validation_data(cooler_2Dselector : cooler.core.RangeSelector2D, 
# #                         promoter_enhancer_dataframe: pd.core.frame.DataFrame, 
# #                         modle_dataframe: pd.core.frame.DataFrame, 
# #                         chrom_name : str, 
# #                         start : int = None, 
# #                         end : int = None,
# #                         count_column : int = 1) -> pd.core.frame.DataFrame:
# #     """Notes promoter/enhancer interaction for a validation dataframe and adds its values to the modle dataframe.

# #     Args:
# #         cooler_2Dselector (cooler.core.RangeSelector2D): a cooler matrix selector
# #         promoter_enhancer_dataframe (pd.core.frame.DataFrame): dataframe with promoter/enhancer regions
# #         modle_dataframe (pd.core.frame.DataFrame): dataframe based on data generated by modle
# #         chrom_name (str): chromosome we want to look at
# #         start (int, optional): start of frame we want to look at. Defaults to None.
# #         end (int, optional): end of frame we want to look at. Defaults to None.

# #     Returns:
# #         pd.core.frame.DataFrame: Modle dataframe with respective validation numbers. 
# #     """
    
# #     if start and end:
# #         print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
# #         valid_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
# #     else:
# #         print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
# #         print("        Note: This may eat a lot memory")
# #         valid_dataframe = cooler_2Dselector.fetch((chrom_name))

# #     resolution = valid_dataframe["start1"][0] - valid_dataframe["end1"][0]

# #     valid_dataframe = note_interactions(cooler_2Dselector, promoter_enhancer_dataframe, chrom_name, start, end, resolution)
    
# #     # Swap columns
# #     modle_count_column = valid_dataframe['modle_count']
# #     valid_count_column = valid_dataframe['valid_count']
# #     valid_dataframe['modle_count'] = valid_count_column
# #     valid_dataframe['valid_count'] = modle_count_column
    

# #     print("Counting up interactions from validation data.")
# #     # Collect all counts within 1000 res dataframe to a 5000 one
# #     for row in modle_dataframe.itertuples():
# #         row_index = row[0]
# #         bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
# #         bin1_start = row[bin1_start_index + 1]
# #         bin1_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_end")
# #         bin1_end = row[bin1_end_index + 1]
        
# #         bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
# #         bin2_start = row[bin2_start_index+ 1]
# #         bin2_end_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_end")
# #         bin2_end = row[bin2_end_index + 1]


# #         start = round_up_and_down(bin1_start,5000)


# #         round_towards_bin1_start = lambda param : round_up_and_down(param,5000)[0] == bin1_start
# #         round_towards_bin2_start = lambda param : round_up_and_down(param,5000)[0] == bin2_start

# #         bin2_in_range = valid_dataframe['bin2_start'].map(round_towards_bin2_start)
# #         filtered_dataframe = valid_dataframe.loc[bin2_in_range]

# #         #print(filtered_dataframe)
# #         bin1_in_range = filtered_dataframe['bin1_start'].map(round_towards_bin1_start)
# #         filtered_dataframe = filtered_dataframe.loc[bin1_in_range]

# #         sum_of_pixel = filtered_dataframe['valid_count'].sum()
        
# #         modle_dataframe.at[row_index, 'valid_count'] = sum_of_pixel


# #     print("Validation data added to dataframe")
# #     return modle_dataframe


# # def create_list_of_chrom(
# #                         cooler_file_path : str, 
# #                         bed_file_path : str, 
# #                         chrom_name : str, 
# #                         start : int = None, 
# #                         end : int = None, 
# #                         cooler_valid_file_path : str = None) -> pd.core.frame.DataFrame:




# #     cooler_object, cooler_2Dselector = read_cool_file(cooler_file_path)

# #     chrom_names : list = cooler_object.chromnames

# #     if chrom_name not in chrom_names:
# #         #TODO Raise error
# #         pass
    
# #     bed_df : pd.core.frame.DataFrame = extract_pls_els_from_bed(bed_file_path)
# #     promoter_enhancer_dataframe : pd.core.frame.DataFrame = filter_type_in_dataframe(bed_df)
    
# #     interactions_dataframe : pd.core.frame.DataFrame = note_interactions(cooler_2Dselector, 
# #                                                                         promoter_enhancer_dataframe, 
# #                                                                         chrom_name, 
# #                                                                         start = start, 
# #                                                                         end = end)
    

# #     if cooler_valid_file_path:
# #         cooler_valid_file_path_5000 = cooler_valid_file_path + "::/resolutions/5000"
# #         cooler_valid_file_path_1000 = cooler_valid_file_path + "::/resolutions/1000"

# #         cool_valid_object_5000 : cooler.api.Cooler = read_mcool_file(cooler_valid_file_path_5000)
# #         valid_2D_selector_5000 = cool_valid_object_5000.matrix(balance=False, as_pixels=True, join=True)

# #         cool_valid_object_1000 : cooler.api.Cooler = read_mcool_file(cooler_valid_file_path_1000)
# #         valid_2D_selector_1000 = cool_valid_object_1000.matrix(balance=False, as_pixels=True, join=True)

# #         interactions_dataframe = note_validation_data(valid_2D_selector_1000, promoter_enhancer_dataframe, interactions_dataframe, chrom_name, start, end)

# #     write_to_file_bedpe(interactions_dataframe, "defaultOutputName.txt")

# #     return interactions_dataframe




# def calculate_promoter_enhancer_bins(promoter_dataframe : pd.DataFrame,
#                                     enhancer_dataframe : pd.DataFrame,
#                                     chrom_name : str,
#                                     chrom_size : str,
#                                     resolution : int,
#                                     out_file_name : str = False,):    

#     print(f"Calculating promoter/enhancer bins in {chrom_name}, res {resolution}")

#     number_of_rows = number_of_columns = math.ceil(chrom_size / resolution)
#     number_of_bins = number_of_rows * number_of_columns

#     matrix_shape = (number_of_columns, number_of_rows)

    
#     # * Filter dataframes to only include entries with the chrom we are looking at currently
#     promoter_dataframe = promoter_dataframe[promoter_dataframe['chrom'] == chrom_name]
#     enhancer_dataframe = enhancer_dataframe[enhancer_dataframe['chrom'] == chrom_name]

#     # * Make sure dataframes are sorted in the correct order
#     promoter_dataframe.sort_values(by=['chromStart','chromEnd'])
#     enhancer_dataframe.sort_values(by=['chromStart','chromEnd'])

    
#     if   number_of_rows >= 4_294_967_295:
#         print("Using np.uint64 as coordinate datatype. This should never happen, because no chromosome is big enough. What are you doing?")
#         coordinate_np_type = np.uint64
#     elif number_of_rows >= 65_535:
#         print("Using np.uint32 as coordinate datatype.")
#         coordinate_np_type = np.uint32
#     else:
#         print("Using np.uint16 as coordinate datatype.")
#         coordinate_np_type = np.uint16

#     # 
#     counter_np_type = np.uint16
#     coordinate_np_type = np.uint16

#     numpy_array_start_size = 10_000

#     coordinate_np = np.empty(shape=(numpy_array_start_size,2), dtype = coordinate_np_type)

#     counter_np = np.empty(shape=numpy_array_start_size, dtype=np.uint)
#     # * Index counter indicating where to add data
#     current_numpy_index = 0

#     for prom_row in promoter_dataframe.itertuples():
#         prom_start = prom_row[2]; 
#         prom_end = prom_row[3]
        
#         # * Round numbers to find exact start and end positions of bins
#         # * We find the start of the first bin, and the end of the second bin
#         # * Calculate how many bins the promoter spans

#         prom_start = int(pelg.round_up_and_down(prom_start, resolution)[0])
#         prom_end = int(pelg.round_up_and_down(prom_end, resolution)[1])
#         number_of_prom_bins = (prom_end - prom_start) // resolution


#         for enh_row in enhancer_dataframe.itertuples():
#             enh_start = enh_row[2]; 
#             enh_end = enh_row[3]

#             # * Round numbers to find exact start and end positions of bins
#             enh_start = int(pelg.round_up_and_down(enh_start, resolution)[0])
#             enh_end = int(pelg.round_up_and_down(enh_end, resolution)[1])
            

#             # * If further apart than 3Mbp, break. This assumes that the enhancer dataframe is sorted by enh_start
#             if abs(enh_start - prom_end) > 3_000_000: 
#                 break

#             # * Calculate how many bins the enhoter spans
#             number_of_enh_bins = (enh_end - enh_start) // resolution

#             # * For each combination of bins, add a pixel
#             for pb in range(number_of_prom_bins):
#                 for eb in range(number_of_enh_bins):
#                     # * Current promoter/enhacer bins we are looking at
#                     prom_bin = prom_start + pb * resolution 
#                     enh_bin = enh_start + eb * resolution

#                     # * Find the indexes of the bins
#                     prom_index = prom_bin // resolution 
#                     enh_index = enh_bin // resolution 

#                     # * Because the matrix is a two dimensional contact matrix, each of the two triangles
#                     # * that make up the matrix (if we split along the diagonal) are flipped duplicates of
#                     # * each other. Therefore we can choose only one of them to keep count of. 
#                     # * If a point appears on the opposite triangle, we instead flip the indexes
#                     # * so all points are on the same triangle. 
#                     indexes = [prom_index, enh_index]
#                     if enh_index < prom_index:
#                         indexes = [enh_index, prom_index]

#                     # * Convert to numpy array before comparison to avoid runtimewarning: invalid value encountered in long_scalars
#                     indexes = np.array(indexes,dtype=coordinate_np_type)
#                     coordinate_index = False
#                     for i, e in enumerate(coordinate_np):
#                         if np.all(e == indexes):
#                             coordinate_index = i
#                             break;

#                     if coordinate_index:
#                         # * If coordinate already exists in array, simply update the counter for the coordinate
#                         counter_np[coordinate_index] = counter_np[coordinate_index] + 1
#                     else:
#                         if current_numpy_index >=  counter_np.size:

#                             new_rows = np.empty((numpy_array_start_size), dtype=counter_np_type)
#                             counter_np = np.append(counter_np, new_rows)


#                             new_rows = np.empty((numpy_array_start_size, coordinate_np.shape[1]), dtype=coordinate_np_type)
#                             coordinate_np = np.vstack((coordinate_np, new_rows))
#                         # * Add values to numpy arrays

#                         coordinate_np[current_numpy_index][0] = indexes[0]
#                         coordinate_np[current_numpy_index][1] = indexes[1]
#                         counter_np[current_numpy_index] = 1
                        
#                         # * Increment index by one 
#                         current_numpy_index += 1

                        
#                         #coordinate_np = np.append(coordinate_np,np.array([indexes]),axis=0)
#                         #counter_np = np.append(counter_np,1)

#     try:
#         # * Resize coordinate_np and counter_np to save memory. New size: current_numpy_index 
#         coordinate_np = coordinate_np[:current_numpy_index]
#         counter_np = counter_np[:current_numpy_index]
#         total_bins_with_regulatory_interaction = len(counter_np)


#         summed_count = np.sum(counter_np,dtype=np.uint)
#         min_count = np.min(counter_np)
#         max_count = np.min(counter_np)
#         median_count = np.median(counter_np)
#         std_count = np.std(counter_np)
#         average_non_zero = summed_count / total_bins_with_regulatory_interaction
#         average_regulatory_count_per_bin_total = summed_count / number_of_bins

#     except RuntimeWarning:
#         print("Runtimewarning in:", chrom_name, resolution )
    
#     dictionary = {'chrom':[chrom_name],'resolution':[resolution],'total_bins_in_chrom':[number_of_bins],
#                 'total_bins_with_pls_and_els':[total_bins_with_regulatory_interaction],
#                 'min_count':[min_count],'max_count':[max_count],'average_count':[average_regulatory_count_per_bin_total],
#                 'median_count':[median_count],'standarddeviation':[std_count],'total_count':[summed_count],
#                 'list_of_counts':[counter_np],'list_of_indexes':[coordinate_np]}
    
#     # * Convert dict to dataframe
#     dataframe = pd.DataFrame.from_dict(dictionary)

#     if not out_file_name:
#         out_file_name = f'./output/predicted_chromosome_data_{chrom_name}_{resolution}.csv'

#     save_dataframe( dataframe=dataframe,
#                     filepath=out_file_name)
    
#     print(f"RETURNING promoter/enhancer bins in {chrom_name}, res {resolution}")
#     return dataframe