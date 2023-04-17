import bioframe as bf # type:ignore
import pandas as pd
import numpy as np
import files
import plot_utils
import cooler
import Constants
import dataframe_functions

def compare(genehancer_input_file : str, cooler_input_file : str, distance_limit : int):

    genehancer_dataframe = files.load_dataframe(genehancer_input_file, allow_pickle=True, numpy_columns=['pixels','assosiated_genes','assosiation_scores'])

    cooler_dataframe = files.load_dataframe(cooler_input_file, allow_pickle=True)

    print(genehancer_dataframe)

    print(cooler_dataframe)
    
    print(distance_limit)

    print(np.sort(np.unique(genehancer_dataframe["promoter_confidence_score"])))

#* Credit: Roberto Rossini
def stolen_code_from_roberto(encode_input, genehancer_input):
    def import_encode_ccres(path) -> pd.DataFrame:
        df = bf.read_table(path, schema="bed").rename(columns={"strand": "feature"})

        df["feature"] = df["feature"].str.replace(r".*PLS.*", "Promoter", regex=True)
        df["feature"] = df["feature"].str.replace(r".*ELS.*", "Enhancer", regex=True)
    
        return df[(df["feature"] == "Promoter") | (df["feature"] == "Enhancer")]

    def import_genehancer_gff(path) -> pd.DataFrame:
        return bf.read_table(path, schema="gff", comment="#")

    df1 = import_encode_ccres(encode_input)
    df2 = import_genehancer_gff(genehancer_input)


    df1 = df1[df1["chrom"].str.contains(r"^chr[\dX]+$", regex=True)]
    df2 = df2[df2["chrom"].str.contains(r"^chr[\dX]+$", regex=True)]

    print(df1)
    print(df2)

    def feature_aware_coverage(df1, df2) -> pd.DataFrame:
        dfs = []
        
        for (feat1, grp1), (feat2, grp2) in zip(sorted(df1.groupby("feature")), sorted(df2.groupby("feature"))):
            assert feat1 == feat2
            dfs.append(bf.coverage(grp1, grp2))
            
        cov = pd.concat(dfs)
        cov["rel_coverage"] = cov["coverage"] / (cov["end"] - cov["start"])

        return cov
        

    def feature_aware_closest(df1, df2) -> pd.DataFrame:
        dfs = []
        
        for (feat1, grp1), (feat2, grp2) in zip(sorted(df1.groupby("feature")), sorted(df2.groupby("feature"))):
            assert feat1 == feat2
            dfs.append(bf.closest(grp1, grp2))
            
        return pd.concat(dfs)

    cov = feature_aware_coverage(df1, df2)
    print(cov)

    df3 = feature_aware_closest(df1, df2)
    print(df3)

    print(df3["distance"].describe().apply(lambda x: format(x, '.2f')))

    


def assosiate_encode_genehancer_id(genehancer_path : str, pls_els_path : str, assosiation_distance : int):

    #print(genehancer_path)

    genehancer_gff_dataframe = pd.read_csv(genehancer_path,delimiter="\t")


    #print(pls_els_path)

    #print(Constants.DATAFRAME_COLUMNS_BED)

    bed_dataframe = pd.read_csv(pls_els_path,delim_whitespace=True,header=None,names=Constants.DATAFRAME_COLUMNS_BED)

    print(len(bed_dataframe))
    print(len(genehancer_gff_dataframe))
    exit()
    #print(len(bf.read_table(pls_els_path, schema="bed").rename(columns={"strand": "type"})))
    #exit()
    pls_dataframe,els_dataframe = dataframe_functions.split_df_to_pls_els(bed_dataframe)

    #print(len(pls_dataframe) + len(els_dataframe))

    # print(pls_dataframe)
    # print(els_dataframe)


    genehancer_gff_promoter_dataframe = genehancer_gff_dataframe[genehancer_gff_dataframe["feature name"] == "Promoter"]
    genehancer_gff_enhancer_dataframe = genehancer_gff_dataframe[genehancer_gff_dataframe["feature name"] == "Enhancer"]
    # genehancer_gff_promoter_dataframe = genehancer_gff_promoter_dataframe.rename(columns={"feature name": "feature"})
    # genehancer_gff_enhancer_dataframe = genehancer_gff_enhancer_dataframe.rename(columns={"feature name": "feature"})

    #print(len(genehancer_gff_promoter_dataframe) + len(genehancer_gff_enhancer_dataframe))

    exit()
    number_of_assosiations_array = np.array([])

    for row in pls_dataframe.itertuples():
        #break # ! remove
        chrom_index = pls_dataframe.columns.get_loc("chrom") + 1
        name_index = pls_dataframe.columns.get_loc("name") + 1
        start_index = pls_dataframe.columns.get_loc("chromStart") + 1
        end_index = pls_dataframe.columns.get_loc("chromEnd") + 1
        
        bed_chrom = row[chrom_index]
        bed_name = row[name_index]
        bed_start = row[start_index]
        bed_end = row[end_index]


        chrom_promoter_dataframe = genehancer_gff_promoter_dataframe[genehancer_gff_promoter_dataframe["#chrom"] == bed_chrom]
        
        min_start = bed_start - assosiation_distance
        max_end = bed_end + assosiation_distance

        associated = chrom_promoter_dataframe[(chrom_promoter_dataframe["end"] >= min_start) & (chrom_promoter_dataframe["start"] <= max_end)]

        if not associated.empty:
            number_of_assosiations_array = np.append(number_of_assosiations_array,len(associated))

    for row in els_dataframe.itertuples():
        #break # ! remove
        chrom_index = els_dataframe.columns.get_loc("chrom") + 1
        name_index = els_dataframe.columns.get_loc("name") + 1
        start_index = els_dataframe.columns.get_loc("chromStart") + 1
        end_index = els_dataframe.columns.get_loc("chromEnd") + 1
        
        bed_chrom = row[chrom_index]
        bed_name = row[name_index]
        bed_start = row[start_index]
        bed_end = row[end_index]


        chrom_enhancer_dataframe = genehancer_gff_enhancer_dataframe[genehancer_gff_enhancer_dataframe["#chrom"] == bed_chrom]
        
        min_start = bed_start - assosiation_distance
        max_end = bed_end + assosiation_distance

        associated = chrom_enhancer_dataframe[(chrom_enhancer_dataframe["end"] >= min_start) & (chrom_enhancer_dataframe["start"] <= max_end)]

        if not associated.empty:
            number_of_assosiations_array = np.append(number_of_assosiations_array,len(associated))
    
    print(np.unique(number_of_assosiations_array))
    print(np.median(number_of_assosiations_array))
    print(np.mean(number_of_assosiations_array))

    print(len(number_of_assosiations_array))

    # def feature_aware_closest(df1, df2) -> pd.DataFrame:
    #     dfs = []
        
    #     for (feat1, grp1), (feat2, grp2) in zip(sorted(df1.groupby("feature")), sorted(df2.groupby("feature"))):
    #         assert feat1 == feat2
    #         dfs.append(bf.closest(grp1, grp2))
        
    #     return pd.concat(dfs)

    # print(pls_dataframe)
    # print(genehancer_gff_promoter_dataframe)

    # df3 = feature_aware_closest(genehancer_gff_promoter_dataframe, pls_dataframe)
    # print(df3)

def aggregation_plot_maybe():
    import seaborn as sns #type:ignore
    import matplotlib.pyplot as plt #type:ignore

    # Load the tips dataset
    tips = sns.load_dataset("tips")

    # Create boxplot with aggregation
    sns.boxplot(x="day", y="total_bill", hue="smoker", data=tips, palette="Set3", dodge=True, linewidth=1, width=0.8, showfliers=False)
    sns.stripplot(x="day", y="total_bill", hue="smoker", data=tips, jitter=True, color=".3", dodge=True, size=4)

    # Add legend and labels
    plt.legend(loc="upper right")
    plt.xlabel("Day")
    plt.ylabel("Total Bill")
    plt.title("Box Aggregation Plot")

    # Show the plot
    plt.savefig("AGGREGATION_TEST.png")
    plt.show()
    plot_utils.show_plot("AGGREGATION_TEST.png")
