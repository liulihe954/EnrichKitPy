## EnrichKitPy

### Introduction

This repo contains off-line python implementation of [EnrichKit Web App](https://github.com/liulihe954/EnrichKitWeb).  
Users are strongly encouraged to read the introductions & instructions from the README of [EnrichKit Web App](https://github.com/liulihe954/EnrichKitWeb).  
Both the web app and the python scripts are internally built on the [EnrichKitDB](https://github.com/liulihe954/EnrichKitDB) database object (available via [Zenodo](https://doi.org/10.5281/zenodo.10535657)).

---

### Usage

- **Parameter file** - `config.py`
  - This file contains necessary parameters for the provided analysis. Users are strongly encouraged to go through the parameters relevant to the desired functions.
  - In order to simply the syntax for each functions, we are gathering all paraters in this single parameter file. So, please note that this composite parameter file may contain parameters that are note relevant in your analysis, you may just ignore them.
- **Function scripts** - `*.py`
  - `id_conversion.py`: This script performs gene ID conversion. To use it, please run `python id_conversion.py input_id_conversion.txt`
  - `loci_annotation.py`: To use it, please run `python loci_annotation.py input_loci_annotation.txt`
  - `pval_aggregation.py`: To use it, please run `python pval_aggregation.py input_pval_aggregation.txt`
  - `enrichment.py`: To use it, please run `python enrichment.py input_enrichment.txt xxxx`; where `xxxx` is the methods to use depending on the nature of the input data. For example, `python enrichment.py input_ora.txt ora` performs Over-Representation Analysis and `python enrichment.py input_prerank.txt prerank` performs GSEA where genes are ranked by `(log2(fc))*-log10(pValue)`.
- **Example input** - `example_input/`
  - In this folder, I am provding input file examples for each of the functions above.
