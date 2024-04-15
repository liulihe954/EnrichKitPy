import pandas as pd
from scipy.stats import combine_pvalues
import statsmodels.stats.multitest as smt
from collections import defaultdict
import scipy.stats as stats
import gseapy as gp

### ### ### ### ### ### ### ### ###
###   for pval aggregation      ###
### ### ### ### ### ### ### ### ###
def pvalue_aggregation(collapse_list, methods):
    
    def create_container():
        container = {"raw_pvals": []}
        for method in methods:
            container[method] = None
        return container
    
    container_dict = defaultdict(create_container)

    for item in collapse_list:
        if not item[0] == "":
            container_dict[item[0]]["raw_pvals"].append(item[1])
    
    for k, v in container_dict.items():
        for method in methods:
            container_dict[k][method] = aggregate(container_dict[k]['raw_pvals'], method)
            
    flattened_data = defaultdict(list)

    for outer_key, inner_dict in container_dict.items():
        flattened_data['Gene'].append(outer_key)
        for inner_key, value in inner_dict.items():
            flattened_data[inner_key].append(value)

    output_df = pd.DataFrame(flattened_data)

    return output_df[['Gene'] + [col for col in output_df.columns if col not in ['Gene', 'raw_pvals']] + ['raw_pvals']]

def aggregate(pvals, method):

    pvals = pd.Series(pvals)

    if method == 'fisher':
        return combine_pvalues(pvals, method='fisher', weights = None)[1]
        
    elif method == 'sidak':
        return min(smt.multipletests(pvals, method = 'sidak')[1])

    elif method == 'simes':
        return simes(pvals)

    elif method == 'fdr':
        return min(smt.multipletests(pvals, method = 'fdr_bh')[1])

    return
    
def simes(p_vals):
    # For a gene containing n p-value, the combined P-value is defined as
    # p{s}=min{np{r}/r;r=1,2…,n} where the p{r} are the individual nucleotide P-values sorted
    # in increasing order. This provides weak control of the family-wise error rate
    # across the set of null hypotheses for all nucleotides in the gene. 
    sorted = p_vals.sort_values()
    ranks = pd.Series(range(1, len(p_vals) + 1))
    multiplied = sorted * len(ranks)
    results = multiplied/ranks.values
    return min(results)







### ### ### ### ### ### ### ### ###
###         for enrichment      ###
### ### ### ### ### ### ### ### ###

def perform_ora(total_genes, sig_genes, candidate_gene_sets, tf = False):

    if not tf:
        # Containers
        db_gene_pool = set()  # Use a set to avoid duplicate entries automatically
        candidate_gene_sets_dict = defaultdict(list)

        # preprocess
        for id, description, gene_member, db_source in candidate_gene_sets:
            if len(gene_member) > 0:
                temp_key = (id, description, db_source)

                if gene_member not in candidate_gene_sets_dict[temp_key]:
                    candidate_gene_sets_dict[temp_key].append(gene_member)
                    db_gene_pool.add(gene_member)

        # fisher's exact test
        N = len(set(total_genes).intersection(db_gene_pool))
        S = len(set(sig_genes).intersection(db_gene_pool))
        DB_Loss_total = (len(total_genes) - N) / len(total_genes)
        DB_Loss_sig = (len(sig_genes) - S) / len(sig_genes)
        output_db = []

        for k, v in candidate_gene_sets_dict.items():
                tmp = ora_helper(k, v, S, N, DB_Loss_total, DB_Loss_sig, total_genes, sig_genes)
                if len(tmp) > 0:
                    output_db.append(tmp)
        return output_db
    
    elif tf:

        # process input gene sets (convert to dictionary)
        # (already done in the main, the ones passed in are converted)

        # total_genes
        # sig_genes


        # process gene dictionary
        db_gene_pool = set()
        candidate_gene_sets_dict = defaultdict(list)

        # preprocess
        for db_source, temp_tf, target_genes, temp_type in candidate_gene_sets:
            if db_source != 'ALL':
                temp_key = (db_source, temp_tf, temp_type)
                temp_gene_list = target_genes.split(',')
                candidate_gene_sets_dict[temp_key] = temp_gene_list
            elif db_source == 'ALL':
                db_gene_pool = target_genes.split(',')

        # perform ora test
        N = len(set(total_genes).intersection(db_gene_pool))
        S = len(set(sig_genes).intersection(db_gene_pool))
        DB_Loss_total = (len(total_genes) - N) / len(total_genes)
        DB_Loss_sig = (len(sig_genes) - S) / len(sig_genes)
        output_db = []

        for k, v in candidate_gene_sets_dict.items():
                tmp = ora_helper(k, v, S, N, DB_Loss_total, DB_Loss_sig, total_genes, sig_genes)
                if len(tmp) > 0:
                    output_db.append(tmp)
        return output_db


def ora_helper(k, v, S, N, DB_Loss_total, DB_Loss_sig, total_gene, sig_gene):

    m = len(set(total_gene).intersection(v))
    if m > 4:
        findG = set(sig_gene).intersection(v)
        s = len(findG)
        cur_table = [[s, S - s], [m - s, N - m - S + s]]
        _, p_value = stats.fisher_exact(cur_table)
        if k[1] and s > 2 and p_value < 1:
            cur_pathway_out = [
                k[0],
                k[2],
                k[1],
                DB_Loss_total,
                DB_Loss_sig,
                s,
                m,
                round(s / m, 2),
                round(p_value, 10),
                ",".join(findG)]
            return cur_pathway_out
    return []

def perform_prerank(rank_list_raw,candidate_gene_sets, tf = False):
    
    if not tf:
        # format gene-sets
        candidate_gene_sets_dict = defaultdict(list)
        for item in candidate_gene_sets:
            temp_key = (item[0],item[1],item[3])
            if len(item[2]) > 1:
                candidate_gene_sets_dict[temp_key].append(item[2])

    elif tf:

        # process gene dictionary
        candidate_gene_sets_dict = defaultdict(list)

        # preprocess
        for db_source, temp_tf, target_genes, temp_type in candidate_gene_sets:
            if db_source != 'ALL':
                temp_key = (db_source, temp_tf, temp_type)
                temp_gene_list = target_genes.split(',')
                candidate_gene_sets_dict[temp_key] = temp_gene_list

    # format gene list
    gene_ranks_df = pd.DataFrame(rank_list_raw)
    gene_ranks_df.set_index(0, inplace=True)

    # run prerank 
    gsea_results_raw = gp.prerank(
        rnk = gene_ranks_df, 
        gene_sets = candidate_gene_sets_dict, 
        min_size = 4,
        permutation_num=2000,
        seed = 1968,
        outdir=None,
        verbose=False
        )
    
    # format output
    gsea_results_out = gsea_results_raw.res2d.sort_values(by=['fdr']).reset_index().rename(columns={'level_0': 'Term','level_1':'Name','level_2':'Source'})
    
    return gsea_results_out

### ### ### ### ### ### ### ### ###
### for processing query results ###
### ### ### ### ### ### ### ### ###

def query_post_process(all_results, curr_function):
    

    if curr_function == 'id_convert':
        all_results = pd.DataFrame(all_results)
        all_results.columns = [
            'query_id',
            'gene_id',
            'ensembl_symbol',
            'entrez_id',
            'ncbi_symbol',
            'vgnc_id',
            'vgnc_symbol',
            'hgnc_orthologs',
            'human_gene_id',
            'human_entrez_id',
            'hgnc_symbol']
        return all_results
    
    elif curr_function == 'enrichment':
        results_other = []
        results_tf = []
        for item in all_results:
            if item[3] == 'tf':
                results_tf.append(item)
            else:
                results_other.append(item)
        return results_other, results_tf
    
    else:
        return all_results

## ### ### ### ### ### ### ### ###
###       miscellaneous         ###
### ### ### ### ### ### ### ### ###

def apply_fdr_correction(group):
    # Perform FDR correction
    _, corrected_pvals, _, _ = smt.multipletests(group['pvalue'], alpha=0.05, method='fdr_bh')
    group['FDR'] = corrected_pvals  # Add corrected p-values to a new column 'FDR'
    return group

def swap_columns(df, col1, col2):
    
    # Ensure the column names are in the DataFrame
    if col1 not in df.columns or col2 not in df.columns:
        raise ValueError("One or both column names do not exist in the DataFrame.")
    
    # Get a list of all column names
    columns = df.columns.tolist()
    
    # Find the index of the two columns to swap
    col1_index, col2_index = columns.index(col1), columns.index(col2)
    
    # Swap the column names in the list
    columns[col1_index], columns[col2_index] = columns[col2_index], columns[col1_index]

    # Use the reordered list of column names to reorder the DataFrame columns
    return df[columns]

def annotateFindG(dictionary_df,target_df):

    # format mapping dict        

    dictionary_df_filtered = dictionary_df[dictionary_df['hgnc_symbol'] != '']
    dictionary_df_filtered.set_index('hgnc_symbol', inplace=True)
    mapping_dict = dictionary_df_filtered['gene_id'].to_dict()

    # modify the target df
    
    target_df['findG'] = target_df.apply(lambda row: match_gene_id(row, mapping_dict), axis=1)

    return target_df

def match_gene_id(row, mapping_dict):
    annotated_list = []
    temp_findG_list = row['findG'].split(',')
    for item in temp_findG_list:
        annotated_list.append(f'{item}={mapping_dict[item]}')    
    return '/'.join(annotated_list)