import argparse
import csv
from random import sample
import config
from src.Setup import Setup
from src.Query import posi_annotate

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
    print('About to process input data.')
    input_tuples = []
    with open(args.input_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            input_tuples.append(tuple(row[0].split(':')))
    print(f"Done. There are {len(input_tuples)} to be annotated.")
    print('')
    
    # for test purposes
    input_list = sample(input_tuples, 400)
    
    # start annotation process
    print('About to annotate the input data.')
    result_df = posi_annotate(config.upstream_length,
                              config.downstream_length,
                              config.splice_length,
                              input_list,
                              config.species)
    print('')
          
    # Output Processing
    print('About to output results.')
    result_df['Input Chr'] = result_df['Input Chr'].astype(int)
    result_df['Input Posi'] = result_df['Input Posi'].astype(int)
    output = result_df.sort_values(by=["Input Chr", "Input Posi"], ascending=[True, True])
    output.drop_duplicates().to_csv("annotation_" + args.input_file, index=False)
    print(f'Done. Results are output to {"annotation_" + args.input_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and annotate genomic positions.")
    parser.add_argument('input_file',
                        help='Input file name for position annotatoin.',
                        type=str)
    args = parser.parse_args()
    main(args)
