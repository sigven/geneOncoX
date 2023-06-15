&nbsp;

# geneOncoX <a href="https://sigven.github.io/geneOncoX/"><img src="man/figures/logo.png" align="right" height="130" width="113"/></a>

__Which human genes are implicated in tumor development__? 

**geneOncoX** is an R package that address this question through the integration of a number of resources with respect to cancer gene annotations, including, but not limited to:

-   [IntOGen](https://www.intogen.org/download) - compendium of mutational cancer driver genes
-   [Network of Cancer Genes](http://ncg.kcl.ac.uk/) - collection of curated cancer genes
-   [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - text-mined predictions of tumor suppressor genes, proto-oncogenes and cancer drivers
-   [Cancer Gene Census](https://cancer.sanger.ac.uk/census) - manually curated resource on cancer genes (soma and germline)
-   [DNA repair genes](https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html) - collection of genes involved in DNA repair
-   [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/) - collections of cancer gene panels used in clinical diagnostics
-   [TSO500 targets](https://emea.illumina.com/products/by-brand/trusight-oncology/tso-500-portfolio.html) - cancer genes targeted by Illumina's TSO500 gene panel
-   [F1CDx targets](https://www.foundationmedicine.com/test/foundationone-cdx) - cancer genes targeted by Foundation One's F1CDx gene panel

The package offers a few pre-processed datasets, along with metadata, that the user can retrieve and use for their own projects or set-ups. The package utilizes the [googledrive](https://googledrive.tidyverse.org/) R package to download the pre-processed and documented datasets to a local cache directory provided by the user.

### Installation

`remotes::install_github('sigven/geneOncoX')`

### Usage

The package offers (currently) five different functions, that each retrieves a specific dataset that can be of use for gene annotation purposes.

-   [`get_basic()`](https://sigven.github.io/geneOncoX/reference/get_basic.html) - retrieves basic, non-transcript-specific gene annotations. Includes tumor suppressor gene/oncogene/driver annotations from multiple resources, NCBI gene summary descriptions, as well as multiple predictions/scores when it comes to gene indispensability and loss-of-function tolerance

-   [`get_gencode()`](https://sigven.github.io/geneOncoX/reference/get_gencode.html) - retrieves two datasets ( *grch37* and *grch38* ) with human gene transcripts from GENCODE, including cross-references to RefSeq, [UniProt](https://www.uniprot.org), [APPRIS](https://appris.bioinfo.cnio.es/#/), and [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)

-   [`get_alias()`](https://sigven.github.io/geneOncoX/reference/get_alias.html) - retrieves a list of gene synonyms, indicating which synonyms are ambiguous or nonambiguous (with respect to primary gene symbols)

-   [`get_predisposition()`](https://sigven.github.io/geneOncoX/reference/get_predisposition.html) - retrieves a list of genes of relevance for cancer predisposition, utilizing multiple resources, including Cancer Gene Census, Genomics England PanelApp, [TCGA's PanCancer study](https://pubmed.ncbi.nlm.nih.gov/29625052/), and others.

-   [`get_panels()`](https://sigven.github.io/geneOncoX/reference/get_panels.html) - retrieves a collection of \> 40 different panels for various cancer conditions, as found in the Genomics England PanelApp.

Technically, each dataset comes as a `list` object in R with

-   a `metadata` data frame that lists URLs, citations, and versions of underlying resources
-   a `records` data frame that contains the actual gene/transcript annotations

### IMPORTANT NOTE

If you use the datasets provided with **geneOncoX**, make sure you properly cite the original publications of the resources integrated, and that you comply with the licensing terms:

1.  IntOGen - [Martínez-Jiménez et al., Nat Rev Cancer, 2020](https://pubmed.ncbi.nlm.nih.gov/32778778/) - [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)
2.  CancerMine - [Lever et al., Nat Methods, 2019](https://pubmed.ncbi.nlm.nih.gov/31110280/) - [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)
3.  Network of Cancer Genes - [Repana et al., Genome Biol, 2019](https://pubmed.ncbi.nlm.nih.gov/30606230/) - *Open Access*
4.  Cancer Gene Census - [Sondka et al., Nat Rev Cancer, 2018](https://pubmed.ncbi.nlm.nih.gov/30293088/) - *Free for non-commercial, academic use - for commercial usage see [https://cancer.sanger.ac.uk/cosmic/license](https://cancer.sanger.ac.uk/cosmic/license)*
5.  DNA repair genes database - [Woods et al., Science, 2001](https://pubmed.ncbi.nlm.nih.gov/11181991/) - *Open Access*
6.  dbNSFP - [Liu et al., Genome Med, 2020](https://pubmed.ncbi.nlm.nih.gov/33261662/) - *Open Access*
7.  Genomics England PanelApp - [Martin et al., Nat Genet, 2019](https://pubmed.ncbi.nlm.nih.gov/31676867/) - *See [licensing terms](https://panelapp.genomicsengland.co.uk/media/files/GEL_-_PanelApp_Terms_of_Use_December_2019.pdf)*
8.  GENCODE - [Frankish et al., Nucleic Acids Res, 2021](https://pubmed.ncbi.nlm.nih.gov/33270111) - *Open Access*

### Contact

sigven AT ifi.uio.no
