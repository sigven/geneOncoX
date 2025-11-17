# Get rank of genes according to association to cancer (Open Targets Platform)

Downloads and returns a dataset that ranks genes according to aggregated
association scores between genes and cancer phenotype terms from the
Open Targets Platform

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a data frame with gene scores/ranks (one record per gene
  and primary site)

## Usage

``` r
get_otp_rank(cache_dir = NA, force_download = FALSE)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

## Value

**metadata** - A data frame with 10 rows and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used

- *abbreviation* - abbreviation used in column names of records

**records** - A data frame with 20,954 rows and 6 columns:

- *entrezgene* - NCBI Entrez gene identifier

- *primary_site* - Primary tumor site

- *tissue_assoc_score* - tissue-specific aggregated association score

- *tissue_assoc_rank* - tissue-specific cancer rank

- *global_assoc_score* - pan-cancer aggregated association score

- *global_assoc_rank* - cancer gene rank (pan-cancer, across tissues)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geneOncoX)
gene_otp_rank <- get_otp_rank(cache_dir = tempdir())
} # }
```
