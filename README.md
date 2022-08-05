## geneOncoX - cancer gene annotations

The goal of **geneOncoX** is to offer an R package that simplifies the process of gene annotation in cancer sequencing projects. The package offers a few pre-processed datasets, along with metadata, that the user can retrieve and use for their own projects or set-ups. Technically, the package utilizes the [googledrive](https://googledrive.tidyverse.org/) R package to download the pre-processed and documented datasets to a local cache directory provided by the user.

**geneOncoX** integrates a number of existing resources when it comes to cancer-relevant gene annotations, including, but not limited to:

-   [IntOGen](https://www.intogen.org/download) - catalog mutational driver genes
-   [Network of Cancer Genes](http://ncg.kcl.ac.uk/) - collection of curated cancer genes
-   [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - text-mined predictions of tumor suppressor genes, proto-oncogenes and cancer drivers
-   [Cancer Gene Census](https://cancer.sanger.ac.uk/census) - manually curated resource on cancer genes (soma and germline)
-   [DNA repair genes](https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html) - collection of genes involved in DNA repair
-   [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/) - collections of cancer gene panels used in clinical diagnostics

### Installation

`remotes::install_github('sigven/geneOncoX')`

### Usage

The package offers currently five different functions, that each retrieves a specific dataset that can be of use for gene annotation purposes.

-   `get_basic()` - retrieves basic, non-transcript-specific gene annotations. Includes tumor suppressor gene/oncogene/driver annotations from multiple resources, NCBI gene summary descriptions, as well as multiple predictions/scores when it comes to gene indispensability and loss-of-function tolerance

-   `get_gencode()` - retrieves two datasets ( *grch37* and *grch38* ) with human gene transcripts from GENCODE, including cross-references to RefSeq, [UniProt](https://www.uniprot.org), [APPRIS](https://appris.bioinfo.cnio.es/#/), and [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)

-   `get_alias()` - retrieves a list of gene synonyms, indicating which synonyms are ambiguous or nonambiguous (with respect to primary gene symbols)

-   `get_predisposition()` - retrieves a list of genes of relevance for cancer predisposition, utilizing multiple resources, including Cancer Gene Census, Genomics England PanelApp, [TCGA's PanCancer study](https://pubmed.ncbi.nlm.nih.gov/29625052/), and others.

-   `get_panels()` - retrieves a collection of \> 40 different panels for various cancer conditions, as found in Genomics England PanelApp

Each dataset comes with a `metadata` object that lists URLs, citations, and versions used.

### IMPORTANT NOTE

If you use the dataset provided with this package, make sure you properly cite the original publications of the resources integrated in **geneOncoX**, i.e.:

1.  IntOGen - [Martínez-Jiménez et al., Nat Rev Cancer, 2020](https://pubmed.ncbi.nlm.nih.gov/32778778/)
2.  CancerMine - [Lever et al., Nat Methods, 2019](https://pubmed.ncbi.nlm.nih.gov/31110280/)
3.  Network of Cancer Genes - [Repana et al., Genome Biol, 2019](https://pubmed.ncbi.nlm.nih.gov/30606230/)
4.  Cancer Gene Census - [Sondka et al., Nat Rev Cancer, 2018](https://pubmed.ncbi.nlm.nih.gov/30293088/)
5.  DNA repair genes database - [Woods et al., Science, 2001](https://pubmed.ncbi.nlm.nih.gov/11181991/)
6.  dbNSFP - [Liu et al., Genome Med, 2020](https://pubmed.ncbi.nlm.nih.gov/33261662/)
7.  Genomics England PanelApp - [Martin et al., Nat Genet, 2019](https://pubmed.ncbi.nlm.nih.gov/31676867/)
8.  GENCODE - [Frankish et al., Nucleic Acids Res, 2021](https://pubmed.ncbi.nlm.nih.gov/33270111)

### Contact

sigven AT ifi.uio.no
