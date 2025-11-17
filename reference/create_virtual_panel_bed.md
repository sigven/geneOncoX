# Function to write virtual panel BED files (CPSR)

Function to write virtual panel BED files (CPSR)

## Usage

``` r
create_virtual_panel_bed(
  gwas_bed_fpath = NA,
  gox_gencode = NULL,
  gox_panels = NULL,
  gox_predisposition = NULL,
  super_panel = FALSE,
  build = "grch37",
  dest_dir = NA
)
```

## Arguments

- gwas_bed_fpath:

  path to assembly-specific GWAS BED track

- gox_gencode:

  object returned by
  [`geneOncoX::get_gencode()`](https://sigven.github.io/geneOncoX/reference/get_gencode.md)

- gox_panels:

  object returned by
  [`geneOncoX::get_panels()`](https://sigven.github.io/geneOncoX/reference/get_panels.md)

- gox_predisposition:

  object returned by
  [`geneOncoX::get_predisposition()`](https://sigven.github.io/geneOncoX/reference/get_predisposition.md)

- super_panel:

  logical indicating if a BED files for a CPSR super panel are to be
  generated

- build:

  genome assembly build (grch37/grch38)

- dest_dir:

  destination directory for BED file

## Value

integer indicating success (0) or failure (-1)
