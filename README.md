# ecDNA in cSCC - Currently Under Construction
This repo contains code for "Extrachromosomal DNA in cutaneous squamous cell carcinoma is associated with increased nodal disease" by 
[Bencomo et al 2024](https://www.biorxiv.org/content/10.1101/2024.02.04.578845v1.abstract). 

## WGS Processing with AmpliconArchitect
A snakemake pipeline was created to run the AmpliconSuite (AmpliconArchitect and associated programs) in parallel on the WGS samples. 
The pipeline can be found in `workflow/`.

Software requirements:
* GATK - the script uses [this](https://hub.docker.com/repository/docker/tbencomo/gatk-bwa-samtools/general) container
* [AmpliconSuite](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) - see github repo for install instructions

Note that the workflow pipeline is memory and CPU intensive - minimum 48G and 16 CPU cores is recommended

Results from each sample can be collected with `gather_results.py` via:

```
gather_results.py [results directory] [output CSV file]
```

## Analysis Scripts

* `get_natgen_genes.py` - extract genes from ecDNA amplicons described in Kim 2020 Nature Genetics paper
* `natgen_lookup.R` - analyze frequency of our ecDNA amplicon genes in Kim 2020 data. Make figures 1a and 1b
* `expression_analysis.R` - compare expression of ecDNA amplicon genes in ecDNA- vs ecDNA+ samples. Make figure 1c
* `create_deseq_data.R` - create DESeq object for `expression_analysis.R`
* `clinical_associations.R` - test clinical features associated with ecDNA status. Make figures 1d-f (Figure 2 in preprint) and code for Table 1 p-values
