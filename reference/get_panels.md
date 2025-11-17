# Get cancer gene panel collections from PanelApp

Downloads and returns a collection of \>40 cancer gene panels from
[Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/).
The dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a data frame with genes found in each panel

## Usage

``` r
get_panels(cache_dir = NA, force_download = FALSE)
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

- *version* - version used

- *abbreviation* - abbreviation used in column names of records

**records** - A data frame with 2,566 rows and 13 columns:

- *genome_build* - human assembly build (grch37/grch38)

- *id* - panel identifier (local)

- *entrezgene* - Entrez gene identifier

- *gene_biotype* - gene biotype

- *genename* - Gene name

- *ensembl_gene_id* - Ensembl gene identifier

- *gepa_moi* - mechanism of inheritance (Genomics England PanelApp)

- *gepa_penetrance* - penetrance (Genomics England PanelApp)

- *gepa_confidence_level* - confidence level (Genomics England PanelApp)

- *gepa_panel_name* - panel name (Genomics England PanelApp)

- *gepa_panel_id* - panel identifier (Genomics England PanelApp)

- *gepa_panel_version* - panel version (Genomics England PanelApp)

- *gepa_phenotype* - associated cancer phenotypes (Genomics England
  PanelApp)

- *gepa_panel_url* - panel URL (Genomics England PanelApp)

## Details

**NOTE:** Gene panel records are provided per genome build (filter on
column `genome_build` to get a build-specific set of panels)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geneOncoX)
gene_panels <- get_panels(cache_dir = tempdir())
} # }
```
