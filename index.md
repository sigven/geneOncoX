# geneOncoX

 

**Which human genes are implicated in tumor development**?

**geneOncoX** is an R package that address this question through the
integration of a number of resources with respect to cancer gene
annotations, including, but not limited to:

- [IntOGen](https://www.intogen.org/download) - compendium of mutational
  cancer driver genes
- [Network of Cancer Genes](http://ncg.kcl.ac.uk/) - collection of
  curated cancer genes
- [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - text-mined
  predictions of tumor suppressor genes, proto-oncogenes and cancer
  drivers
- [CanVar-UK](https://canvaruk.org/) - cancer predisposition genes
- [DNA repair
  genes](https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html) -
  collection of genes involved in DNA repair
- [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/) -
  collections of cancer gene panels used in clinical diagnostics
- [TSO500
  targets](https://emea.illumina.com/products/by-brand/trusight-oncology/tso-500-portfolio.html) -
  cancer genes targeted by Illumina’s TSO500 gene panel
- [F1CDx
  targets](https://www.foundationmedicine.com/test/foundationone-cdx) -
  cancer genes targeted by Foundation One’s F1CDx gene panel

The package offers a few pre-processed datasets, along with metadata,
that the user can retrieve and use for their own projects or set-ups.
The package utilizes the
[googledrive](https://googledrive.tidyverse.org/) R package to download
the pre-processed and documented datasets to a local cache directory
provided by the user.

### Getting started

- [Installation
  instructions](https://sigven.github.io/geneOncoX/articles/geneOncoX.html#installation)
- [Usage
  examples](https://sigven.github.io/geneOncoX/articles/geneOncoX.html#get-basic-gene-annotations)

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/sigven/geneOncoX/blob/main/.github/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
