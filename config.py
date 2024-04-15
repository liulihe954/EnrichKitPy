import os

### Environment Variables ###
BASE = os.getcwd()
cloud_store = 'zenodo'
zenodo_get = 'zenodo_get'
zenodo_record_id = '10535657'

### Analyses Parameters ###

## Species in the current analysis
# please use lower case abbreviation 
species = "bta" # cap, ovi, gal, sus, equ

# Other supported species
# 'bta' -> 'Cow - Bos Taurus - ARS-UCD1.2'
# 'cap' -> 'Goat - Capra_hircus - ARS1'
# 'ovi' -> 'Sheep - Ovis_aries - Oar_v3.1'
# 'gal' -> 'Chicken - Gallus_gallus - gca000002315v5.GRCg6a'
# 'sus' -> 'Pig - Sus_scrofa - Sscrofa11.1'
# 'equ' -> 'Horse - Equus_caballus - EquCab3.0'

## Thr provided gene ID type for the ID conversion function
input_id_type = 'ensembl' # or 'entrez', current supports these two


## Genomic coordinates parameters for loci annotation & pvalue aggregationx
upstream_length = 10000
downstream_length = 10000
splice_length = 50 # to define the length of splice donor and accepter
aggregation_methods = ['fisher','sidak','simes','fdr']


## Parameters for enrichment analysis - ORA & GSEA
enrichment_input_id_type = 'ensembl' ## Note: currently only accept 'ensembl' gene id, the above "input_id_type" will be ignored.
pathway_database = ['tf','kegg','go','interpro','mesh','reactome']

## All supported databases, please use lower case abbreviation 
# 'go' -> 'Gene Ontology(GO)'
# 'interpro' -> 'Interpro'
# 'kegg' -> 'KEGG'
# 'mesh'  -> 'Medical Subject Headings (MeSH)'
# 'reactome' -> 'Reactome'
# 'tf' -> 'Human Transcription Factor'                                                            




