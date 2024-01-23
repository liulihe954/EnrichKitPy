import os
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import sqlite3
import pandas as pd
from .QueryExecutor import QueryExecutor
from .SQLs import SQLs

########
######## posi_annotate()
########

def posi_annotate(input_upstream, input_downstream, input_splice, input_query_list, input_species):
    queries = input_query_list
    all_results = []

    # Thread-safe progress bar initialization
    with tqdm(total=len(queries), desc='Processing Queries') as pbar:
        # Determine the number of cores to use
        max_cores = max(1, os.cpu_count() - 2)

        # Using ThreadPoolExecutor for parallel queries with max_cores
        with ThreadPoolExecutor(max_workers=max_cores) as thread_executor:
            # Map the process_query function to all queries
            futures = {thread_executor.submit(process_query, query, input_upstream, input_downstream, input_splice, input_species): query for query in queries}
            
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    all_results.extend(result)  # Collect results
                pbar.update(1)  # Update the progress bar

    # Creating the result DataFrame
    result_df = pd.DataFrame(
        all_results,
        columns=['Input Chr', 'Input Posi', 'Label', 'Info 1', 'Info 2', 'Gene ID', 'Gene Biotype', 'Strand'])

    return result_df


########
######## process_query()
########

def process_query(cur_gene_params_raw, input_upstream, input_downstream, input_splice, input_species):
    
    SPLICE_OFFSET = input_splice - 50
    cur_gene_params =  tuple([input_species]) + cur_gene_params_raw
    
    all_results = []
    
    try:
        
        # Create a new executor for each query
        executor = QueryExecutor(
            '.dbcache/EnrichKitDB.sqlite',
            upstream=input_upstream,
            downstream=input_downstream
        )
        
        ##
        ultra_limit = False

        # s1: query target gene
        cur_gene_results = executor.query(SQLs['posi2gene'], category = 'posi2gene', params=cur_gene_params)

        if cur_gene_results:
            for item in cur_gene_results:

                gene_target = executor.query(SQLs['ekid2geneid'], category = 'ekid2geneid', params=(str(item[0]),))

                ############################ With GENE TARGET ################################################

                if gene_target:

                    cur_posi_params = (
                        gene_target[0][0], # current ek_gene_id
                        int(cur_gene_params[2]), # current posi
                    )

                    ###########################################################################
                    # altra upstream OR downstream
                    ###########################################################################

                    # sure to return something this line
                    gene_target_start, gene_target_end, gene_target_strand, gene_target_gene_biotype = executor.query(SQLs['getGenelimit'], category = 'getGenelimit', params=(str(gene_target[0][0]),))[0]

                    if cur_posi_params[1] <= (gene_target_start-10000):
                        ultra_limit = True
                        ultra_limit_formt = [
                            cur_gene_params[1],      # chromosome
                            cur_gene_params[2],      # position
                            'feature',              # label
                            'upstream',
                            '',
                            gene_target[0][2],       # Gene ID
                            gene_target_gene_biotype,# biotype
                            gene_target_strand       # strand
                        ]
                        all_results.append(list(ultra_limit_formt))

                    elif cur_posi_params[1] >= (gene_target_end + 10000):
                        ultra_limit = True
                        ultra_limit_formt = [
                            cur_gene_params[1],      # chromosome
                            cur_gene_params[2],      # position
                            'feature',              # label
                            'downstream',
                            '',
                            gene_target[0][2],       # Gene ID
                            gene_target_gene_biotype,# biotype
                            gene_target_strand       # strand
                        ]
                        all_results.append(list(ultra_limit_formt))

                    if not ultra_limit:
                        ###########################################################################
                        # exon
                        ###########################################################################
                        exon_results = executor.query(SQLs['posi2exon'], category = 'posi2exon', params = cur_posi_params)
                        if exon_results:
                            for exon_item in exon_results:
                                exon_item_format = [
                                    cur_gene_params[1], # chromosome
                                    cur_gene_params[2], # position
                                    'exonic',           # label
                                    exon_item[1],       # Exon id
                                    exon_item[6],       # Transcript ID
                                    gene_target[0][2],  # Gene ID
                                    exon_item[4],       # biotype
                                    exon_item[5]        # strand
                                ]
                                all_results.append(list(exon_item_format))

                        else:
                            pass

                        ###########################################################################
                        # features
                        ###########################################################################
                        feature_results = executor.query(SQLs['posi2feature'], category = 'posi2feature', params = cur_posi_params)
                        if feature_results:

                            for feature_item in feature_results:

                                if feature_item[1] == 'transcript':
                                    feature_item_format = [
                                        cur_gene_params[1], # chromosome
                                        cur_gene_params[2], # position
                                        'feature',          # label
                                        feature_item[1],    # feature name
                                        str(feature_item[2]) + "-" + str(feature_item[3]),# None (should be Transcript ID) 
                                        gene_target[0][2],  # Gene ID
                                        feature_item[4],       # biotype
                                        feature_item[5]        # strand
                                        ]                                
                                else:
                                    feature_item_format = [
                                        cur_gene_params[1], # chromosome
                                        cur_gene_params[2], # position
                                        'feature',          # label
                                        feature_item[1],    # feature name
                                        '',                 # None (should be Transcript ID) 
                                        gene_target[0][2],  # Gene ID
                                        feature_item[4],       # biotype
                                        feature_item[5]        # strand
                                        ]

                                all_results.append(list(feature_item_format))
                        else:              
                            pass

                        ###########################################################################
                        # features 1  "posi2cfeature1"
                        ###########################################################################
                        cfeature_results1 = executor.query(SQLs['posi2cfeature1'], category = 'posi2cfeature1', params = cur_posi_params)

                        if cfeature_results1:
                            for feature_item in cfeature_results1:

                                cfeature_item_format = [
                                    cur_gene_params[1], # chromosome
                                    cur_gene_params[2], # position
                                    'feature',          # label
                                    feature_item[1],    # feature name
                                    '',                 # None (should be Transcript ID) 
                                    gene_target[0][2],  # Gene ID
                                    feature_item[4],       # biotype
                                    feature_item[5]        # strand
                                ]

                                if feature_item[1] == 'splice donor' and int(cur_gene_params[2]) > feature_item[3] + SPLICE_OFFSET:
                                    cfeature_item_format[3] = 'intron'

                                elif feature_item[1] == 'splice acceptor' and int(cur_gene_params[2]) < feature_item[2] - SPLICE_OFFSET:
                                    cfeature_item_format[3] = 'intron'

                                all_results.append(list(cfeature_item_format))

        else:
            loci_annotation_placeholder = [
                cur_gene_params[1], # chromosome
                cur_gene_params[2], # position
                'intergenic or out of bound',
                '',
                '',
                '',
                '',
                '']
            all_results.append(list(loci_annotation_placeholder))
            
        executor.close()    
        return all_results
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    

