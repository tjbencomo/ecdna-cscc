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

