# gwas2gtca
Converting GWAS summary statistics to GCTA COJO analysis input file



### Installation

To install all the dependencies needed in a conda environment, use:

```bash
conda create -n gwas2gcta -c conda-forge -c bioconda r-base r-dplyr r-data.table r-purrr bioconductor-genomicranges bioconductor-rtracklayer
```

Then, once all the dependencies are installed, run

```bash
conda activate gwas2gcta
```



### Example

```bash
Rscript gwas2gcta.R example_gwas.txt CHROM POS_b37 HG37 ALT REF POOLED_ALT_AF EFFECT_SIZE SE N test_output
```

