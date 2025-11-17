# Get human cancer predisposition genes

Downloads and retrieves a pre-processed dataset of human cancer
predisposition genes (CPGs). Genes deemed relevant for cancer
predisposition have been collected from multiple resources:

- [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/)

- [TCGA's PanCancer study](https://pubmed.ncbi.nlm.nih.gov/29625052/)

- [CanVar-UK](https://canvaruk.org/)

- Other/manually curated

The dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a list with records of cancer predisposing gene, one
  record per gene

## Usage

``` r
get_predisposition(cache_dir = NA, force_download = FALSE)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten

## Value

**metadata** - A data frame with 1 row and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used

- *abbreviation* - abbreviation used in column names of records

**records** - A data frame with 594 rows and 8 columns:

- *symbol* - official gene symbol

- *entrezgene* - Entrez gene identifier

- *gene_biotype* - gene biotype

- *cpg_moi* - mechanism of inheritance (AD, AR)

- *cpg_syndrome_cui* - Concept unique identifiers (CUI, UMLS) -
  inherited cancer syndromes

- *cpg_cancer_cui* - Concept unique identifiers (CUI, UMLS) - inherited
  cancer conditions

- *cpg_source* - Sources supporting predisposition gene PANEL_APP,
  CANVAR-UK, TCGA_PANCAN_2018, OTHER, ACMG_SF

- *cpg_phenotypes* - associated cancer phenotypes

- *cpg_mod* - mechanism of disease)

## Details

**NOTE**: The dataset also contains genes recommended for reporting of
incidental findings (ACMG_SF).

## Examples

``` r
if (FALSE) { # \dontrun{
library(geneOncoX)
gene_predisp <- get_predisposition(cache_dir = tempdir())
} # }
```
