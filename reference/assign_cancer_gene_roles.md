# Function to classify genes as tumor suppressors/proto-oncogenes, driver etc. based on multiple lines of evidence.

Function to classify genes as tumor suppressors/proto-oncogenes, driver
etc. based on multiple lines of evidence.

## Usage

``` r
assign_cancer_gene_roles(
  gox_basic = NULL,
  min_citation_support = 5,
  min_citation_support_cm_only = 20,
  min_sources_driver = 2
)
```

## Arguments

- gox_basic:

  output from geneOncoX::get_basic()

- min_citation_support:

  minimum citation count for support from CancerMine

- min_citation_support_cm_only:

  minimum citation count for genes annotated with cancer gene roles from
  CancerMine only

- min_sources_driver:

  minimum number of sources (NCG, CanVar-UK, CancerMine, TCGA, IntoGen)
  that must contribute to cancer driver gene status

## Value

data frame with cancer gene annotation status
