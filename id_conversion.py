import argparse
import csv
import pandas as pd
from random import sample
import config
from src.Setup import Setup
from src.Query import process_query
from src.Utils import pvalue_aggregation

def main(args):
    # Set up Environment
    print('About to check database availability.')
    setup = Setup(config.BASE,
                  config.zenodo_get,
                  config.zenodo_record_id,
                  config.cloud_store)
    setup.download()
    print('')

    # Process Input
    print('Start processing input data.')
    input_list = []

    with open(args.input_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
    #             print(row)
                input_list.append(
                    (row[0],)
                )
    input_list = list(set(input_list))

    print(f"Done. ID type {config.input_id_type}")
    print(f"There are {len(input_list)} gene IDs to be converted.")
    print('')
    
    # for test purposes
    input_list = sample(input_list, 50) # [('5','118503702') for _ in range(1000)] 
    
    # start annotation process
    print('Start ID conversion.')
    input_params = (
    'id_convert',
    config.species,
    config.input_id_type,
    input_list,
    len(input_list)
    )
    result_df = process_query(input_params)
    print('')
          
    # Output Processing
    print('Start to output results.')
    output_file_path = "OUTPUT_ID_conversion_" + args.input_file
    result_df.drop_duplicates().to_csv(output_file_path, index=False)

    print(f'Done. Results are output to {output_file_path}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gene ID conversion.")
    parser.add_argument('input_file',
                        help='Input file contains single column of all the gene IDs (ensembl or entrez) to be converted.',
                        type=str)
    args = parser.parse_args()
    main(args)
