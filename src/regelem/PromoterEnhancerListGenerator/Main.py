from .Functions import *

import Constants
import sys
import time

start_time = 0

def main():

    start_time : float = time.time()

    if len(sys.argv) < OBLIGATORY_ARGUMENT_AMOUNT:
        print("Expected", OBLIGATORY_ARGUMENT_AMOUNT, "arguments. Was given", len(sys.argv))
        exit()

    cooler_file_path = sys.argv[1]
    pls_els_bed_path = sys.argv[2]
    extrusion_barrier_bed_path = sys.argv[3]

    cooler_object, cooler_2Dselector = read_cool_file(cooler_file_path)

    chrom_names = cooler_object.chromnames

    bed_dataframe = extract_pls_els_from_bed(pls_els_bed_path)

    pls_els_dataframe = filter_type_in_dataframe(bed_dataframe)

    note_interactions(cooler_2Dselector, pls_els_dataframe, chrom_names)

    return

def exit_program():
    end_time = time.time()
    runtime = end_time - start_time

    minutes = round(runtime / 60)
    seconds = runtime % 60

    print("Finished process in", minutes, "m", seconds, "s")
    exit()

if __name__ == '__main__':
    main()