import argparse
import csv
import pandas as pd
from random import sample
import config
from src.Setup import Setup
from src.Query import process_query

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
    input_tuples = []
    
    with open(args.input_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            input_tuples.append(tuple(row[0].split(':'))) # TODO: add error handling
                
    print(f"Done. There are {len(input_tuples)} rows to be annotated.")
    print('')
    
    # for test purposes
    # input_list = sample(input_tuples, 50) # [('5','118503702') for _ in range(1000)] 
    
    # start annotation process
    print('Start annotating the input data.')
    input_params = (
        'posi_annotate',
        config.upstream_length,
        config.downstream_length,
        config.splice_length,
        input_list,
        config.species,
        len(input_list)
    )
    all_results = process_query(input_params)
    
    # Creating the result DataFrame
    result_df = pd.DataFrame(
        all_results,
        columns=['Input Chr', 'Input Posi', 'Label', 'Info 1', 'Info 2', 'Gene ID', 'Gene Biotype', 'Strand'])

    print('')
          
    # Output Processing
    print('Start to output results.')
    result_df['Input Chr'] = result_df['Input Chr'].astype(int)
    result_df['Input Posi'] = result_df['Input Posi'].astype(int)
    output = result_df.sort_values(by=["Input Chr", "Input Posi"], ascending=[True, True])
    output.drop_duplicates().to_csv("OUTPUT_annotation_" + args.input_file, index=False)
    print(f'Done. Results are output to OUTPUT_{"annotation_" + args.input_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and annotate genomic positions.")
    parser.add_argument('input_file',
                        help='Input file name for position annotatoin.',
                        type=str)
    args = parser.parse_args()
    main(args)
