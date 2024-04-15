import os
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import sqlite3
import pandas as pd
from .QueryExecutor import QueryExecutor
from .SQLs import SQLs
from .SQLs import sqlite_object
from .Utils import query_post_process


########  
######## process_query layer ########
########
        
def process_query(input_params):
    # initiate constants
    curr_function = input_params[0]
    all_results = []
    
    # Create a single SQLite connection in read-only mode
    conn = sqlite3.connect(sqlite_object, uri=True, check_same_thread=False)
    
    # Thread-safe progress bar initialization, the last element in the input_params is the number of queries to execute
    with tqdm(total=input_params[-1], desc='Processing Queries') as pbar:
        
        # Determine the number of cores to use
        max_cores = max(1, os.cpu_count() - 2)

        # Using ThreadPoolExecutor for parallel queries with max_cores
        with ThreadPoolExecutor(max_workers=max_cores) as thread_executor:    
            # Map the process_query function to all queries
            
            try:
                # for posi_annotate feature
                if curr_function == 'posi_annotate':
                    futures = {thread_executor.submit(process_query_posi_annotate, 
                                                      conn,
                                                      query,
                                                      input_params[1],
                                                      input_params[2],
                                                      input_params[3],
                                                      input_params[5]): query for query in input_params[4]}
                    
                elif curr_function == 'pval_aggregate':
                    futures = {thread_executor.submit(process_query_pval_aggregate, 
                                                      conn,
                                                      query,
                                                      input_params[1],
                                                      input_params[2],
                                                      input_params[4]): query for query in input_params[3]}
                elif curr_function == 'id_convert':
                    futures = {thread_executor.submit(id_convert, 
                                                      conn,
                                                      query,
                                                      input_params[1],
                                                      input_params[2]): query for query in input_params[3]}                
                elif curr_function == 'enrichment':
                    futures = {thread_executor.submit(enrichment, 
                                                      conn,
                                                      query,
                                                      input_params[1]): query for query in input_params[2]}
                    
            except Exception as e:
                print(f'Error in ThreadPoolExecutor: {e}')

                
            # collect results
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    all_results.extend(result)  # Collect results
                pbar.update(1)  # Update the progress bar
            
    conn.close()

    return query_post_process(all_results, curr_function)


######## individual process_query customization layer ########

########
######## enrichment() ########
########
def enrichment(conn, cur_gene_params_raw, input_species):

    cur_gene_params = (input_species, cur_gene_params_raw)

    try:
        executor = QueryExecutor(conn, upstream=0, downstream=0)

        if cur_gene_params_raw == 'tf':

            cur_gene_results = executor.query(SQLs['extract_tf_gene'], category = 'extract_tf_gene', params = cur_gene_params)
        
        else:

            cur_gene_results = executor.query(SQLs['extract_geneset_ensembl'], category = 'extract_geneset', params = cur_gene_params)

        return [sublist + (cur_gene_params_raw,) for sublist in cur_gene_results]
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    
    return None
    
########
######## pval_annotate() ########
########
def process_query_posi_annotate(conn, cur_gene_params_raw, input_upstream, input_downstream, input_splice, input_species):
    
    SPLICE_OFFSET = input_splice - 50
    cur_gene_params =  tuple([input_species]) + cur_gene_params_raw
    
    all_results = []
    
    try:
        
#         # Create a new executor for each query
#         executor = QueryExecutor(
#             '.dbcache/EnrichKitDB.sqlite',
#             upstream=input_upstream,
#             downstream=input_downstream
#         )
        executor = QueryExecutor(conn, upstream=input_upstream, downstream=input_downstream)
        ##
        ultra_limit = False
        
        #
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
            
        return all_results
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    

########
######## id_convert() ########
########
def id_convert(conn, cur_gene_params_raw, input_species, input_id_type):
    #
    cur_gene_params =  (input_species,) + tuple(cur_gene_params_raw) #tuple([input_species]) + cur_gene_params_raw 

    all_results = []

    try:
        executor = QueryExecutor(conn, upstream=0, downstream=0)

        if input_id_type == 'ensembl':
            cur_gene_results = executor.query(SQLs['id_convert_ensembl'], category = 'id_convert', params = cur_gene_params)
            
        elif input_id_type == 'entrez':
            cur_gene_results = executor.query(SQLs['id_convert_entrez'], category = 'id_convert', params = cur_gene_params)
        
        if cur_gene_results:
            temp_output = [cur_gene_params[1],]
            temp_output.extend(cur_gene_results[0][2:])
            all_results.append(temp_output)

        else:
            place_holder = [
                cur_gene_params[1],
            ]
            place_holder.extend(['NA'] * 10)

            all_results.append(place_holder)

        return all_results
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    


########
######## pval_aggregate() ########
########
def process_query_pval_aggregate(conn, cur_gene_params_raw, input_upstream, input_downstream, input_species):
    #
    cur_gene_params =  tuple([input_species]) + cur_gene_params_raw
    all_results = []

    try:
        executor = QueryExecutor(conn, upstream=input_upstream, downstream=input_downstream)

        cur_gene_results = executor.query(SQLs['posi2gene'], category = 'posi2gene', params=cur_gene_params)

        if cur_gene_results:
            for item in cur_gene_results:
                gene_target = executor.query(SQLs['ekid2geneid'], category = 'ekid2geneid', params=(str(item[0]),))
                temp_out = [gene_target[0][2],cur_gene_params[3]]
                all_results.append(temp_out)
        else:
            pval_aggregation_placeholder = [
                '',
                cur_gene_params[3]
                ]
            all_results.append(pval_aggregation_placeholder)
            
        return all_results
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    

    ########
