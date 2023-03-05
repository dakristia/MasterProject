"""Constants """
# FILE_PATH_IN = './data/input/'
# FILE_MODLE_COOL = FILE_PATH_IN + 'output5.cool'
# FILE_PLS_ELS_BED = FILE_PATH_IN + 'PLS_pdELS.7group.bed'
# FILE_EXTRUSION_BARRIER = FILE_PATH_IN + 'GRCh38_GM12878_barriers_RAD21_occupancy.bed'
# FILE_GENE_ANNOTATION = FILE_PATH_IN + 'hg38.refGene.gtf'

FILE_PATH_OUT = './data/output/output.bedpe'

OBLIGATORY_ARGUMENT_AMOUNT = 4


# Index values for BED
BED_CHROM = 0
BED_CHROM_START = 1
BED_CHROM_END = 2
BED_NAME = 3
BED_SCORE = 4
BED_STRAND = 5
BED_THICKSTART = 6
BED_THICKEND = 7
BED_ITEMRGB = 8
BED_BLOCKCOUNT = 9
BED_BLOCKSIZES = 10
BED_BLOCKSTARTS = 11

# BEDPE constants
BEDPE_CHROM = 0
BEDPE_CHROM_START1 = 1
BEDPE_CHROM_END1 = 2
BEDPE_CHROM_START2 = 3
BEDPE_CHROM_END2 = 4
BEDPE_NAME = 5
BEDPE_SCORE = 6
BEDPE_STRAND1 = 7
BEDPE_STRAND2 = 8

# Index values for GTF
GTF_SEQNAME = 0
GTF_SOURCE = 1
GTF_FEATURE = 2
GTF_START = 3
GTF_END = 4
GTF_SCORE = 5
GTF_STRAND = 6
GTF_FRAME = 7
GTF_ATTRIBUTE = 8

#
PROMOTER = "promoter"
ENHANCER = "enhancer"
PLS = "PLS"
ELS = "ELS"



#The dataframe format we'll use within the list generator
DATAFRAME_COLUMNS_INTERNAL = ["chrom", "enh_start", "enh_end", "prom_start", "prom_end", 
                            "enh_name", "prom_name", "bin1_start", "bin1_end", 
                            "bin2_start", "bin2_end", "modle_count", "valid_count"]
DATAFRAME_COLUMNS_COOLER = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]
DATAFRAME_COLUMNS_BED = ["chrom","chromStart","chromEnd","name","score","strand"
                        ,"thickStart","thickEnd","itemRgb","type","complete"]
#DATAFRAME_COLUMNS_BED = ["chrom1", "start1",
#                         "end1", "chrom2", "start2", "end2", "count"]
DATAFRAME_COLUMNS_PLSELS = ["chrom", "p_name", "p_start",
                            "p_end", "e_name", "e_start", "e_end", "count", "freq"]
DATAFRAME_COLUMNS_PLSELS_OW = [
    "chrom", "regType", "regionname"] + DATAFRAME_COLUMNS_BED
DATAFRAME_COLUMNS_BEDPE = ["chrom1", "start1", "end1", "chrom2",
                            "start2", "end2", "name", "score", "strand1", "strand2"]
DATAFRAME_COLUMNS_STATISTICAL = ['distance','mincount','maxcount','averagecount','mediancount',
                                'standarddeviation','totalcount','numberofcounts','allcounts']


DISTANCE_LIMIT = 3_000_000 #3 Mega bp
RESOLUTION_BP = 5000
