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
                    (row[0].split(':')[0],
                    row[0].split(':')[1],
                    float(row[1])
                    )
                )
    input_list = list(set(input_list))
    print(f"Done. There are {len(input_list)} rows to be aggregated.")
    print('')
    
    # for test purposes
    # input_list = sample(input_list, 50) # [('5','118503702') for _ in range(1000)] 
    
    # start annotation process
    print('Start annotating the input data.')
    input_params = (
        'pval_aggregate',
        config.upstream_length,
        config.downstream_length,
        input_list,
        config.species,
        len(input_list)
    )
    collapse_list = process_query(input_params)
    result_df = pvalue_aggregation(collapse_list, config.aggregation_methods)
    print('')
          
    # Output Processing
    print('Start to output results.')
    output_file_path = "OUTPUT_aggregation_" + args.input_file
    result_df['raw_pvals'] = result_df['raw_pvals'].apply(tuple)
    result_df.drop_duplicates().to_csv(output_file_path, index=False)

    print(f'Done. Results are output to {output_file_path}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate pvalues of genomic positions to genes.")
    parser.add_argument('input_file',
                        help='Input file name for pvalue aggregation.',
                        type=str)
    args = parser.parse_args()
    main(args)
