# # # #Pixel operations
   
# # #     # pixelRange = modle_cool.pixels(join=True)
# # #     # print("Pixeldata: ")
# # #     # print(pixelRange.shape)
# # #     # print(type(pixelRange))
# # #     # print(pixelRange[0:30])



# # # #Getting data from dict

# # #     # print(list(PLS_ELS_dict.keys()))
# # #     # for l in PLS_ELS_dict["chr1"][PLS][0:4]:
# # #     #     print(l)
# # #     # print()
# # #     # for l in PLS_ELS_dict["chr1"][ELS][0:4]:
# # #     #     print(l)

# # #     #print(len(PLS_ELS_dict["chr1"][PLS])) 


# # # #Some cooler / matrix operations

# # #     # for col in example_matrix.columns: print(col)
# # #     # example_matrix = modle_cool.matrix(balance=False, as_pixels=True, join=True)[0:5, 0:5]
# # #     # print(example_matrix.columns)


    
# # #     #print("start", start, "end", end, range(start,end))#,  end="\r")
# # #     #print("range")
    
# # #     # print(type(pls_start),type(start))
# # #     # print(type(df))
# # #     # for y in range(start,end):
# # #     #     print(y)

# # #     #TODO: Serves as a note for now. Remove this when you feel like its not needed.
# # #     #mat = c.matrix(balance=False).fetch('chr5')


# # #     #print(c.chromnames)
# # #     # print(modle_cool.info)
# # #     # print(modle_cool.binsize)
# # #     # print(modle_cool.chromnames)
# # #     # print(modle_cool.chromsizes)
# # #     # print(modle_cool.bins())
# # #     # print(modle_cool.pixels(join=True))
# # #     # print(modle_cool.matrix(balance=False, as_pixels=True, join=True)[1000:1005, 1000:1005])
# # #     # print()
# # #     #print(modle_cool.matrix(balance=False, as_pixels=True, join=True)[0:5, 0:5])

# # #     #print(modle_cool.matrix(balance=False, as_pixels=True, join=True)[0:5, 0:5]["start1"][6])

    
    

# # # def extract_PLS_ELS_From_Bed(bedList):
# # #     """Extract a list of all lines containing either PLS, dELS or pELS

# # #     Keyword arguments:
# # #     bedList -- a list containing strings, where each element is a line in a .bed file
# # #             -- a filename string to a .bed file

# # #     Returns: 
# # #     Dict containing all .bed data on PLS and ELS. Dict key is the chromosome name, value is a dict with two keys either 'PLS' or 'ELS'.
# # #     Each value in these two entries is a list of bed data. Same indexes as columns in a bed data. 
# # #     """

# # #     if type(bedList) is str:
# # #         bedList = read_Bed(bedList)



# # #     bed_dict = {}

# # #     for bedLine in bedList:
# # #         chromName = bedLine[BED_CHROM]
# # #         if chromName not in bed_dict.keys():
# # #             bed_dict[chromName] = {PLS: [], ELS: []}


# # #         blockcount = bedLine[BED_BLOCKCOUNT]
# # #         blockcount_split = blockcount.split(",")

# # #         if PLS in blockcount_split: 
# # #             bed_dict[chromName][PLS].append(bedLine)

        
# # #         for word in blockcount_split:
# # #             if re.search(ELS,word):
# # #                 bed_dict[chromName][ELS].append(bedLine)
# # #                 break

# # #     return bed_dict


# # # def extract_PLS_ELS_From_Bed_List(bedList: list) -> dict:
# # #     """Extract a list of all lines containing either PLS, dELS or pELS

# # #     Keyword arguments:
# # #     bedList -- a list containing strings, where each element is a line in a .bed file

# # #     Returns: 
# # #     Dict containing all .bed data on PLS and ELS. Dict key is the chromosome name, value is a dict with two keys either 'PLS' or 'ELS'.
# # #     Each value in these two entries is a list of bed data. Same indexes as columns in a bed data. 
# # #     """

# # #     bed_dict: dict = {}

# # #     for bedLine in bedList:
# # #         chrom_name = bedLine[BED_CHROM]
# # #         if chrom_name not in bed_dict.keys():
# # #             bed_dict[chrom_name] = {PLS: {}, ELS: {}}

# # #         blockcount = bedLine[BED_BLOCKCOUNT]
# # #         blockcount_split = blockcount.split(",")

# # #         if PLS in blockcount_split:
# # #             bed_dict[chrom_name][PLS][bedLine[BED_NAME]] = bedLine

# # #         for word in blockcount_split:
# # #             if re.search(ELS, word):
# # #                 bed_dict[chrom_name][ELS][bedLine[BED_NAME]] = bedLine
# # #                 break

# # #     return bed_dict

# # # def write_to_file_bedpe(dataframe: pd.DataFrame, out_file_path: str) -> None:
# # #     """
# # #     Writes the contents of a dataframe to a .bedpe file. 

# # #     Arguments:
# # #         dataframe: a pandas dataframe. Requires the columns: 'chrom', 'p_name', 'e_name' , 'count', 'freq'
# # #         out_file_path: the file path of the output file
# # #     Returns:
# # #         None

# # #     """

# # #     contact_list = []

# # #     for contact_index, contact_line in dataframe.iterrows():
# # #         contact_list.append([contact_line['chrom'],
# # #                              contact_line['p_start'], contact_line['p_end'],
# # #                              contact_line['e_start'], contact_line['e_end'],
# # #                              contact_line['p_name'] + "-" +
# # #                              contact_line['e_name'],
# # #                              contact_line['freq']])

# # #     print("Writing output to file: ", out_file_path)
# # #     with open(out_file_path, 'w') as out_file:
# # #         for line in contact_list:
# # #             for e in line:
# # #                 out_file.write(str(e) + " ")
# # #             out_file.write("\n")

# # #     return

# # # def get_extent_of_PLSELS(dict, chrom_name):
# # #     #Finding lowest/highest values of PLS and ELS. Get that range from df ---------

# # #     #! Doesn't work as dict format has changed
# # #     lowest = int(dict[chrom_name][PLS][0][BED_CHROM_START])
# # #     highest = lowest
# # #     for row in dict[chrom_name][PLS].values():
# # #         num = int(row[BED_CHROM_START])
# # #         num2 = int(row[BED_CHROM_END])
# # #         if num < lowest:
# # #             lowest = num
# # #         if num2 > highest:
# # #             highest = num

# # #     for row in dict[chrom_name][ELS].values():
# # #         num = int(row[BED_CHROM_START])
# # #         num2 = int(row[BED_CHROM_END])
# # #         if num < lowest:
# # #             lowest = num
# # #         if num2 > highest:
# # #             highest = num

# # #     return (lowest,highest)

# # # #Returns a tuple. First element is a dataframe with the regions, second element is a dictionary with range rounded down to nearest resolution as key
# # # def find_pls_els_regions(dataframe, pls_els_dict, chrom_name):

# # #     df = dataframe
# # #     df_reduced = pd.DataFrame(columns=df.columns)
# # #     df_reduced_pls_in_1 = pd.DataFrame(columns=df.columns)
# # #     df_reduced_pls_in_2 = pd.DataFrame(columns=df.columns)
    
# # #     dict_elsrange_aliases = defaultdict(list)

# # #     # Get all rows that may contain a PLS as either first or second range
# # #     # Add them to two seperate dataframes
# # #     for pls_row in pls_els_dict[chrom_name][PLS].values():
# # #         pls_start = int(pls_row[BED_CHROM_START])
# # #         pls_end = int(pls_row[BED_CHROM_END])

# # #         pls_start_rounded_down = round_up_and_down(pls_start,RESOLUTION_BP)[0]
# # #         pls_start_rounded_up = round_up_and_down(pls_start,RESOLUTION_BP)[1]


# # #         pls_end_rounded_down = round_up_and_down(pls_end,RESOLUTION_BP)[0]

# # #         frame1 = frame2 = frame3 = frame4 = None
# # #         if pls_start_rounded_down == pls_end_rounded_down:
# # #             frame1 = df.loc[df['start1'] == pls_start_rounded_down]
# # #             frame3 = df.loc[df['start2'] == pls_start_rounded_down]
# # #         else:
# # #             frame1 = df.loc[df['start1'] == pls_start_rounded_down]
# # #             frame2 = df.loc[df['start1'] == pls_start_rounded_up]
# # #             frame3 = df.loc[df['start2'] == pls_start_rounded_down]
# # #             frame4 = df.loc[df['start2'] == pls_start_rounded_up]


# # #         df_reduced_pls_in_1 = pd.concat([df_reduced_pls_in_1,frame1,frame2])
# # #         df_reduced_pls_in_2 = pd.concat([df_reduced_pls_in_2,frame3,frame4])
 
# # #     # Get all rows that may contain an ELS as either first or second range
# # #     # Looking for these values in dataframes made by last for-loop
# # #     for els_row in pls_els_dict[chrom_name][ELS].values():
# # #         els_start = int(els_row[BED_CHROM_START])
# # #         els_end = int(els_row[BED_CHROM_END])

# # #         # Round to nearest 5000 to 
# # #         els_start_rounded_down = round_up_and_down(els_start,RESOLUTION_BP)[0]
# # #         els_start_rounded_up = round_up_and_down(els_start,RESOLUTION_BP)[1]

# # #         els_end_rounded_down = round_up_and_down(els_end,RESOLUTION_BP)[0]

# # #         #Temporary dataframes
# # #         frame1 = frame2 = frame3 = frame4 = None

# # #         # If entire els is within same region (based on resolution)
# # #         if els_start_rounded_down == els_end_rounded_down:
# # #             frame1 = df_reduced_pls_in_1.loc[df_reduced_pls_in_1['start2'] == els_start_rounded_down]
# # #             frame3 = df_reduced_pls_in_2.loc[df_reduced_pls_in_2['start1'] == els_start_rounded_down]

# # #             dict_elsrange_aliases[els_start_rounded_down].append(els_row[BED_NAME])
# # #         #If els extends over two regions, we look at both
# # #         else:
# # #             frame1 = df_reduced_pls_in_1.loc[df_reduced_pls_in_1['start2'] == els_start_rounded_down]
# # #             frame2 = df_reduced_pls_in_1.loc[df_reduced_pls_in_1['start2'] == els_start_rounded_up]
# # #             frame3 = df_reduced_pls_in_2.loc[df_reduced_pls_in_2['start1'] == els_start_rounded_down]
# # #             frame4 = df_reduced_pls_in_2.loc[df_reduced_pls_in_2['start1'] == els_start_rounded_up]

# # #             dict_elsrange_aliases[els_end_rounded_down].append(els_row[BED_NAME])

# # #         #Add all data to the same dataframe
# # #         df_reduced = pd.concat([df_reduced,frame1,frame2,frame3,frame4])

# # #     #Remove all duplicates, sort by the first area then the second, reset indexes to match new data
# # #     df_reduced = df_reduced.drop_duplicates().sort_values(by=['start2']).sort_values(by=['start1']).reset_index(drop=True)
 
# # #     return (df_reduced, dict_elsrange_aliases)

# # # def note_interactions(matrix_selector, pls_els_dict, chromnames, small_test = True):
# # #     #TODO: Create a dict with every interaction every chromosome has had
# # #     #TODO: First find out how many interactions this will contain

# # #     #TODO Run this on every chromname
# # #     for chrom_name in chromnames:

# # #         #Entire matrix
# # #         df = matrix_selector.fetch((chrom_name))
        

# # #         n_of_rows = df.shape[0]
# # #         print("rows in df: ",n_of_rows)
# # #         print("row in PLS dict:", len(pls_els_dict[chrom_name][PLS]))
# # #         print("row in ELS dict:", len(pls_els_dict[chrom_name][ELS]))

# # #         # Find lowest and highest places where PLS or ELS are present, and cut away everything outside of this 
# # #         # print(pls_els_dict)
# # #         # lowest, highest = get_extent_of_PLSELS(pls_els_dict,chrom_name,matrix_selector)
# # #         # df = matrix_selector.fetch((chrom_name, lowest, highest))
# # #         # print(df)


# # #         df_reduced, dict_elsrange_aliases = find_pls_els_regions(df,pls_els_dict,chrom_name)

# # #         # It is quite possible that you have a dataframe (df_reduced) with all PLS and ELS interactions
# # #         #TODO: Make a list of these interactions by connecting their names to the values, and adding counts


# # #         df_contact = annotate_promoter_enhancer_interactions(df_reduced, DATAFRAME_COLUMNS_PLSELS, pls_els_dict, dict_elsrange_aliases, chrom_name)

# # #         if small_test: break #If this is on, will only do first chromosome

# # #     return df_contact


# # # def annotate_promoter_enhancer_interactions(target_dataframe, new_dataframe_columns, pls_els_dict, dict_elsrange_aliases, chrom_name):

# # #     df_contact = pd.DataFrame(columns=new_dataframe_columns)

# # #     els_dict = pls_els_dict[chrom_name][ELS]

# # #     total_count = 0
# # #     for count in target_dataframe['count']: total_count += count

# # #     for pls_key, pls_row in pls_els_dict[chrom_name][PLS].items():
# # #         pls_start = int(pls_row[BED_CHROM_START])
# # #         pls_end = int(pls_row[BED_CHROM_END])

# # #         pls_data = [pls_row[BED_CHROM],pls_row[BED_NAME],pls_row[BED_CHROM_START],pls_row[BED_CHROM_END]]

# # #         pls_start_rounded_down = round_up_and_down(pls_start,RESOLUTION_BP)[0]
# # #         pls_end_rounded_down = round_up_and_down(pls_end,RESOLUTION_BP)[0]

# # #         frame1 = target_dataframe.loc[target_dataframe['start1'] == pls_start_rounded_down]
# # #         frame2 = pd.DataFrame()
# # #         if pls_start_rounded_down != pls_end_rounded_down: frame2 = target_dataframe.loc[target_dataframe['start1'] == pls_end_rounded_down]

# # #         df_pls = pd.concat([frame1,frame2])

# # #         df_temp = pd.DataFrame(columns=DATAFRAME_COLUMNS_PLSELS)


# # #         for pd_index, pd_row in df_pls.iterrows():
# # #             els_range = pd_row['start2']

# # #             potential_references = dict_elsrange_aliases[els_range]

# # #             for els_name in potential_references:

# # #                 els_row = els_dict[els_name]
# # #                 els_data = [els_row[BED_NAME],els_row[BED_CHROM_START],els_row[BED_CHROM_END]]

# # #                 column_data = pls_data + els_data + [pd_row['count'], round(((pd_row['count'])/total_count),10)]


# # #                 df_temp = pd.DataFrame([column_data],columns=DATAFRAME_COLUMNS_PLSELS)
# # #                 df_contact = pd.concat([df_contact,df_temp])
        

# # #     df_contact = df_contact.drop_duplicates().sort_values(by=['e_name']).sort_values(by=['p_name']).reset_index(drop=True)


# # #     return df_contact

# # # def main():

# # #     start_time : float = time.time()

# # #     if len(sys.argv) < OBLIGATORY_ARGUMENT_AMOUNT:
# # #         print("Expected", OBLIGATORY_ARGUMENT_AMOUNT, "arguments. Was given", len(sys.argv))
# # #         exit()



# # #     cooler_file_path = sys.argv[1]
# # #     pls_els_bed_path = sys.argv[2]
# # #     extrusion_barrier_bed_path = sys.argv[3]
# # #     print(cooler_file_path)

# # #     modle_cool, cooler_2Dselector = read_cooler_file(cooler_file_path)

# # #     # Getting the .bed data for promoters/enhancers in 
    
# # #     pls_els_bed_list = read_Bedfile(pls_els_bed_path)

# # #     #Not used?
# # #     extrusion_barrier_list = read_Bedfile(extrusion_barrier_bed_path)

# # #     print("Extracting PLS_ELS")
# # #     pls_els_dict = extract_PLS_ELS_From_Bed_List(pls_els_bed_list)
    

# # #     print("Noting interactions")
# # #     df_contact = note_interactions(cooler_2Dselector, pls_els_dict, modle_cool.chromnames)

# # #     write_to_file_bedpe(df_contact,FILE_PATH_OUT)

# # #     exit_program()

# # #     #TODO: Consider if this should be kept

# # #     print("Extrusion barrier:") 

# # #     for line in extrusion_barrier_list[:30]: print(line)


# # #     print()
# # #     gtf_list = read_Gtfile(FILE_GENE_ANNOTATION)

# # #     print("GTF: ")
# # #     for line in gtf_list[:20]: print(line)



# # #     extrusion_barrier_dict = extract_Extrusion_Barrier_From_Bed(extrusion_barrier_list)


# # # def note_validation_data(cooler_2Dselector, modle_dataframe: pd.core.frame.DataFrame, chrom_name : str, start : int = None, end : int = None) -> pd.core.frame.DataFrame:
    
    
# # #     if start and end:
# # #         print("Fetching data from dataframe related to", chrom_name, "starting at", start, "ending at", end)
# # #         valid_dataframe = cooler_2Dselector.fetch((chrom_name,start,end))
# # #     else:
# # #         print("No start and/or end specified. Fetching entirety of chromosome: ", chrom_name)
# # #         print("        Note: This may eat a lot memory")
# # #         valid_dataframe = cooler_2Dselector.fetch((chrom_name))


# # #     print("\n\n\n\n\n DATAFRAME")
# # #     print(modle_dataframe)
# # #     print(valid_dataframe)
# # #     print("\n\n\n\n\n DATAFRAME")

# # #     #print(modle_dataframe)#.loc[valid_dataframe["start1"] == modle_dataframe["bin1_start"]])
# # #     #print(valid_dataframe)

# # #     #df.query('C1>=0 and C1<=20 and C2>=0 and C2<=20 and C3>=0 and C3<=20')
# # #     #print(valid_dataframe.query('start1 == 60000 and start2 == 65000'))

# # #     #df[(df>=0)&(df<=20)].dropna()

# # #     mixed_dataframe = pd.core.frame.DataFrame()

# # #     for row in modle_dataframe.itertuples():
# # #         bin1_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin1_start")
# # #         bin1_start = row[bin1_start_index]
# # #         bin2_start_index = DATAFRAME_COLUMNS_INTERNAL.index("bin2_start")
# # #         bin2_start = row[bin2_start_index]

# # #         #valid_dataframe.loc[valid_dataframe["start1"] == bin1_start]

# # #         renamed_frame = modle_dataframe.rename(columns={"bin1_start":"start1", "bin2_start":"start2"})
# # #         #print(renamed_frame)
# # #         mergedStuff = pd.merge(valid_dataframe, renamed_frame, on=['start1','start2'], how='inner')
# # #         #print(mergedStuff)

# # #         valid_count_series = mergedStuff['count'].rename('valid_count')

# # #         #print(valid_count_series[0])
# # #         modle_dataframe["valid_count"] = valid_count_series
# # #         modle_dataframe = modle_dataframe.fillna(-1).astype({"valid_count":"int"})

# # #     print(modle_dataframe)
# # #     return modle_dataframe


# # # * FINISHED
# # def read_extrution_barrier_file(chrom_name="chr19", 
# #                                 resolution=1000,
# #                                 start=None,
# #                                 end=None,
# #                                 extrution_barrier_bed_file="chr19_extrution_barriers.bed",
# #                                 chrom_sizes_file = "H1.chrom.sizes"):
# #     extrution_barrier_bed_file_path = input_folder + extrution_barrier_bed_file
# #     chrom_size_file_path = input_folder + chrom_sizes_file

# #     extrution_barrier_df : pd.core.frame.DataFrame = pd.read_csv(extrution_barrier_bed_file_path, 
# #                                                 delim_whitespace=True, 
# #                                                 index_col=False,
# #                                                 header = None, 
# #                                                 names = ["chrom","chromStart","chromEnd","name","score","strand"])



# #     chrom_size_dataframe : pd.core.frame.DataFrame = pd.read_csv(chrom_size_file_path,
# #                                                     delim_whitespace=True,
# #                                                     index_col="chrom",
# #                                                     header = None,
# #                                                     names = ["chrom", "size"])


# #     extrution_barrier_dict = extrution_barrier_df.to_dict()#.pop('chrom').pop('name').pop('score')
# #     del extrution_barrier_dict['chrom']; del extrution_barrier_dict['name']; del extrution_barrier_dict['score']

# #     chrom_size = chrom_size_dataframe.at['chr19','size']


# #     rows = chrom_size // resolution + 1

# #     array = np.zeros((rows,))
# #     array = np.full(rows,"",dtype=str)

    

# #     for i in range(0,len(extrution_barrier_dict['strand'])):
# #         chromStart = extrution_barrier_dict['chromStart'][i]
# #         chromEnd = extrution_barrier_dict['chromEnd'][i]
# #         strand = extrution_barrier_dict['strand'][i]

# #         chrom_start_bin = chromStart // resolution
# #         chrom_end_bin = chromEnd // resolution

# #         if strand not in array[chrom_start_bin]: array[chrom_start_bin] += strand
# #         if strand not in array[chrom_end_bin]: array[chrom_end_bin] += strand
        
# #     array = np.char.replace(array, '+-','b')
# #     array = np.char.replace(array,'-+','b')

# #     start_index = start // resolution
# #     end_index = end // resolution

# #     array = array[start_index:end_index+1]

# #     import matplotlib.pyplot as plt
# #     import sys
# #     import subprocess

# #     img = np.zeros((1, len(array), 3), dtype=np.uint8)

# #     # Set the color of each pixel based on the value of the corresponding array element
# #     img[0,np.where(array == '')] = [128, 128, 128]
# #     img[0,np.where(array == '+')] = [255, 0, 0]
# #     img[0,np.where(array == '-')] = [0, 0, 255]

# #     # Plot the image
# #     fig, ax = plt.subplots()
# #     ax.imshow(img, aspect=10, extent=(0, len(array), 0, 1))
# #     ax.set_xticks([])
# #     ax.set_yticks([])
# #     ax.spines['top'].set_visible(False)
# #     ax.spines['right'].set_visible(False)
# #     ax.spines['bottom'].set_visible(False)
# #     ax.spines['left'].set_visible(False)


# #     plt.savefig('./testlineplot.png')

# #     imageViewerFromCommandLine = {'linux':'xdg-open',
# #                                     'win32':'explorer',
# #                                     'darwin':'open'}[sys.platform]
# #     subprocess.run([imageViewerFromCommandLine, './testlineplot.png'])

# # def corr_coef():

# #     matrix1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
# #     matrix2 = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])

# #     correlation_coefficient = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0,1]
# #     print("Correlation coefficient:", correlation_coefficient)


# def generate_plotFigure_for_TADs_section_of_essay():
#     #bed_file_new = test_folder + "2022-11-07_cisRegElements_ENCFF573NKX_ENCFF760NUN_ENCFF919FBG_ENCFF332TNJ_grepped.7group.bed"
#     hic_data_file_path = input_folder + "outfile_binsize5000.cool"
#     valid_data_file_path = input_folder + "4DNFI9GMP2J8.mcool::/resolutions/5000"
#     #hic_data_file_path = input_folder + "OUTPUT.cool"

#     chrom_name = "chr19"
#     start = "400_000"
#     end = "1_000_000"
#     #start = None
#     #end = None
#     #start = "0"
#     #end = "1_000_000"
#     #dataframe = pelg.Dataframe_Functions.extract_original_dataframe_from_cool(hic_data_file_path, chrom_name, start, end)

#     cooler_object : cooler.api.Cooler = cooler.Cooler(hic_data_file_path)
#     valid_cooler_object : cooler.api.Cooler = cooler.Cooler(valid_data_file_path)
    
#     #send dataframe to plot function.
#     coolerplotter = plotter.CoolerPlotter()

#     if start and end:
#         modle_annotation="MoDLE " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
#         valid_annotation="H1 " + chrom_name + ":" + start.replace("_",",") + "-" + end.replace("_",",")
#     else:
#         modle_annotation="MoDLE " + chrom_name
#         valid_annotation="H1 " + chrom_name

#     coolerplotter.simple_contact_matrix_plot(in_data = cooler_object,chrom_name = chrom_name,start=start,end=end,title="MoDLE counts",annotation=modle_annotation,out_filepath="figure_TAD_section_essay.png",open_in_viewer=True)
#     coolerplotter.simple_contact_matrix_plot(in_data = valid_cooler_object,chrom_name = chrom_name,start=start,end=end,title="valid counts",annotation=valid_annotation,out_filepath="figure_valid_TAD_section_essay.png",open_in_viewer=True)
# #plot_distance_line_plot("chr19")
# # 
