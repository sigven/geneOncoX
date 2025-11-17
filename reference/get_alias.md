# Get human gene aliases from NCBI

Downloads and returns a dataset that indicate ambiguous and unambiguous
gene aliases/synonyms for human genes. Gene aliases with less than three
characters have been ignored, and a few custom aliases have been added
(source = `custom`). The dataset comes as a `list` object, with two
elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a list with gene aliases indicating ambiguous/
  non-ambiguous state

## Usage

``` r
get_alias(cache_dir = NA, force_download = FALSE)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

## Value

**metadata** - A data frame with 1 row and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used/datestamp

- *abbreviation* - abbreviation used in column names of records

**records** - A data frame with 177,916 rows and 5 columns:

- *alias* - gene alias/synonym

- *symbol* - primary symbol

- *entrezgene* - Entrez gene identifier

- *n_primary_map* - number of primary symbols linked to the alias

- *ambiguous* - logical indicating if alias is ambiguous or not

- *source* - source for gene synonyms (NCBI, custom)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geneOncoX)
gene_alias <- get_alias(cache_dir = tempdir())
} # }
```
