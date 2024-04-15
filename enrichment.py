import argparse
import csv
import pandas as pd
from random import sample
import config
from src.Setup import Setup
from src.Query import process_query
from collections import defaultdict
from src.Utils import perform_ora, apply_fdr_correction, swap_columns, perform_prerank, annotateFindG
import numpy as np


def main(args):

    # Set up Environment
    print('About to check database availability.')
    setup = Setup(config.BASE,
                    config.zenodo_get,
                    config.zenodo_record_id,
                    config.cloud_store)
    setup.download()
    print('')

    ### Process input
    enrichment_type  = args.type
    input_file = args.input_file

    if enrichment_type == 'ora':
        total_genes = []
        sig_genes = []
        with open(input_file, 'r') as file:
                reader = csv.reader(file)
                for row in reader:
                    total_genes.append(row[0])
                    if row[1] == '1':
                        sig_genes.append(row[0])
        if 'tf' in config.pathway_database:
            print('Internal Gene ID conversion will firstly be performed for human TF enrichment')
            ensembl_hgnc_mapping = process_query((
                'id_convert',
                config.species,
                config.enrichment_input_id_type,
                [(item,) for item in total_genes],
                len(total_genes)))[['gene_id','hgnc_symbol']]
            
            total_genes_hgnc = ensembl_hgnc_mapping['hgnc_symbol'].to_list()

            total_genes_hgnc_filter = list(filter(None, total_genes_hgnc))

            sig_genes_hgnc = process_query((
                'id_convert',
                config.species,
                config.enrichment_input_id_type,
                [(item,) for item in sig_genes],
                len(sig_genes)))['hgnc_symbol'].tolist()
            
            sig_genes_hgnc_filter = list(filter(None, sig_genes_hgnc))

    elif enrichment_type == 'prerank':
        rank_list_raw = []
        all_gene = set()
        with open(input_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                rank_list_raw.append(
                    [row[0],float(row[1]) * -np.log10(float(row[2]))]
                )
                all_gene.add(row[0])

        if 'tf' in config.pathway_database:
            print('Internal Gene ID conversion will firstly be performed for human TF enrichment')
            # container 
            rank_list_raw_hgnc = []
            # obtain a mapping
            ensembl_hgnc_mapping = process_query((
                'id_convert',
                config.species,
                config.enrichment_input_id_type,
                [(item,) for item in all_gene],
                len(all_gene)))[['gene_id','hgnc_symbol']]
            
            dictionary_df_filtered = ensembl_hgnc_mapping[ensembl_hgnc_mapping['gene_id'] != '']
            dictionary_df_filtered.set_index('gene_id', inplace=True)
            mapping_dict = dictionary_df_filtered['hgnc_symbol'].to_dict()

            # rank_list_raw_hgnc
            for item in rank_list_raw:
                temp_alternative = mapping_dict[item[0]]
                if temp_alternative != '':
                    hgnc_alternative = [temp_alternative,item[1]]
                    rank_list_raw_hgnc.append(hgnc_alternative)
            
    else:
        print('Please correct the enrichment type parameter, currently only ora or prerank is available')


    ### Perform enrichment
    print('Start retrieving pathway list from the database.')
    input_params = (
        'enrichment',
        config.species,
        config.pathway_database,
        len(config.pathway_database)
    )

    # one query to take out all the relevant pathways (no need to branch for ora/gsea, multi-threads)
    candidate_gene_sets, candidate_gene_sets_tf = process_query(input_params)


    print(f'Start performing enrichment analysis using {enrichment_type} method for databases:{config.pathway_database}')
    # run enrichment
    if enrichment_type == 'ora':
        if len(candidate_gene_sets) != 0:
            output_df = pd.DataFrame(
                perform_ora(total_genes, sig_genes, candidate_gene_sets),
                columns=[
                    'Term',
                    'Source',
                    'Name',
                    'DB_Loss_total',
                    'DB_Loss_sig',
                    'Num sigG',
                    'Num totalG',
                    'hitsPerc',
                    'pvalue',
                    'findG'
                ])
            ora_output_fdr = output_df.groupby('Source', as_index=False).apply(apply_fdr_correction).reset_index(drop=True).sort_values(by=['FDR','pvalue'])
            ora_output_fdr_swapped = swap_columns(ora_output_fdr, 'findG', 'FDR')

            ora_output_fdr_swapped.drop_duplicates().to_csv(f"OUTPUT_{enrichment_type}_enrichment_{args.input_file}", index=False)
            print(f"Done. Results for {[item for item in config.pathway_database if item != 'tf']} are output to OUTPUT_{enrichment_type}_enrichment_{args.input_file}")


        if len(candidate_gene_sets_tf) != 0:
            print(f'Now performing enrichment for TF using {enrichment_type}')
            ora_output_df_tf = pd.DataFrame(
                perform_ora(total_genes_hgnc_filter, sig_genes_hgnc_filter, candidate_gene_sets_tf, tf = True),
                columns=[
                    'Source',
                    'placeholder',
                    'TF_Name',
                    'DB_Loss_total',
                    'DB_Loss_sig',
                    'Num sigG',
                    'Num totalG',
                    'hitsPerc',
                    'pvalue',
                    'findG'
                ])
            ora_output_df_tf_fdr = ora_output_df_tf.groupby('Source', as_index=False).apply(apply_fdr_correction).reset_index(drop=True).sort_values(by=['FDR','pvalue'])
            ora_output_df_tf_fdr_swapped = swap_columns(ora_output_df_tf_fdr, 'findG', 'FDR')
            ora_output_df_tf_fdr_swapped_annotateFindG = annotateFindG(ensembl_hgnc_mapping, ora_output_df_tf_fdr_swapped)

            ora_output_df_tf_fdr_swapped_annotateFindG.drop_duplicates().to_csv(f"OUTPUT_{enrichment_type}_TF_enrichment_{args.input_file}", index=False)
            print(f'Done. TF enrichment results are output to OUTPUT_{enrichment_type}_TF_enrichment_{args.input_file}')

    elif enrichment_type == 'prerank':
        if len(candidate_gene_sets) != 0:
            output_df = perform_prerank(rank_list_raw,candidate_gene_sets)
            
            output_df.drop_duplicates().to_csv(f"OUTPUT_{enrichment_type}_enrichment_{args.input_file}", index=False)
            print(f"Done. Results for {[item for item in config.pathway_database if item != 'tf']} are output to OUTPUT_{enrichment_type}_enrichment_{args.input_file}")


        if len(candidate_gene_sets_tf) != 0:
            print(f'Now performing enrichment for TF using {enrichment_type}')
            output_df_tf = perform_prerank(rank_list_raw_hgnc, candidate_gene_sets_tf, tf = True)
            output_df_tf.rename(columns={'genes': 'findG'}, inplace=True)
            output_df_tf_annotateFindG = annotateFindG(ensembl_hgnc_mapping, output_df_tf)
            
            output_df_tf_annotateFindG.drop_duplicates().to_csv(f"OUTPUT_{enrichment_type}_TF_enrichment_{args.input_file}", index=False)
            print(f'Done. TF enrichment results are output to OUTPUT_{enrichment_type}_TF_enrichment_{args.input_file}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perfrom gene set enrichment with ORA or prerank.")
    parser.add_argument('input_file',
                        help='Input file name for position annotation.',
                        type=str),
    parser.add_argument('type',
                        choices=["ora", "prerank"],
                        help='The type of enrichment analysis, ora or prerank.',
                        type=str),
    args = parser.parse_args()
    main(args)