import pandas as pd
import numpy as np
import cooler


import sys
sys.path.append('../src/regelem')
import Constants

def main():

    cooler_object1 = cooler.Cooler("../input/H1-hESC.mcool::/resolutions/5000")
    cooler_object2 = cooler.Cooler("../input/H1-hESC.mcool::/resolutions/1000")

    selector1 = cooler_object1.matrix(balance=False, as_pixels=True, join=True)
    selector2 = cooler_object2.matrix(balance=False, as_pixels=True, join=True)  

    offset = 10_000

    df1 = selector1.fetch(("chr1",1_000_000 + offset,1_010_000 + offset))
    df2 = selector2.fetch(("chr1",1_000_000 + offset,1_010_000 + offset))

    print(df1.query(f'start1 < 1_015_000 and start2 >= 1_015_000'))
    print(df2.query(f'start1 < 1_015_000 and start2 >= 1_015_000'))

    print(np.sum(df2.query(f'start1 < 1_015_000 and start2 >= 1_015_000')['count']))

if __name__ == "__main__":

    main()


def old_test():
    dataframe = pd.read_csv("../input/H1-hESC.7group.PLSELSONLY.bed",delim_whitespace=True,header=None,names=Constants.DATAFRAME_COLUMNS_BED)

    #dataframe = dataframe[dataframe["chrom"] == "chr1"]
    pls_dataframe = dataframe[dataframe["type"].str.contains('PLS')]
    els_dataframe = dataframe[dataframe["type"].str.contains('ELS')]

    chrom_names = np.unique(np.array(dataframe["chrom"]))

    count = 0

    for name in chrom_names:
        chrom_dataframe = dataframe[dataframe["chrom"] == name]
        pls_dataframe = dataframe[dataframe["type"].str.contains('PLS')]
        els_dataframe = dataframe[dataframe["type"].str.contains('ELS')]

        # start_column = np.array(chrom_dataframe["chromStart"])
        # end_column = np.array(chrom_dataframe["chromEnd"])
        print("number of promoters:", len(pls_dataframe))
        print("number of enhancers:", len(els_dataframe))

        print()
        for e in pls_dataframe.itertuples():
            for y in els_dataframe.itertuples():
                if abs(e[2] - y[3]) > 2_999_900: break
                if abs(e[3] - y[2]) > 2_999_900: break
                if abs(e[2] - y[3]) < 50 or abs(e[3] - y[2]) < 50 : count += 1
                #count += 1
                # if abs(e - y) <= 100:
                #     print(y, e)

        print(f"count after {name}: {count}")
        exit()
    print(count)
    exit()

    sizes = np.array([])

    for row in chrom_dataframe.itertuples():
        chromStart = row[2]
        chromEnd = row[3]
        size = abs( chromStart - chromEnd)
        #if chromStart % 100 == 50 and chromEnd % 100 == 50: print(row)
        #if size == 150: print(row)
        sizes = np.append(sizes,size )

    print(np.min(sizes))
    print(np.max(sizes))
    print(np.mean(sizes))
    print(np.median(sizes))

    print(np.count_nonzero())

    print(len(sizes))
