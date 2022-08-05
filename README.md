## geneOncoXREF - cancer gene annotations

The goal of **geneOncoXREF** is to offer an R package that simplifies the process of gene annotation in cancer sequencing projects. The package offers a few pre-processed datasets, along with metadata, that the user can retrieve and use for their own projects or set-ups. Technically, the package utilizes the [googledrive](https://googledrive.tidyverse.org/) R package to download the pre-processed and documented datasets to a local cache directory provided by the user.

**geneOncoXREF** integrates a number of existing resources when it comes to cancer-relevant gene annotations, including, but not limited to:

-   [IntOGen](https://www.intogen.org/download) - catalog mutational driver genes
-   [Network of Cancer Genes](http://ncg.kcl.ac.uk/) - collection of curated cancer genes
-   [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - text-mined predictions of tumor suppressor genes, proto-oncogenes and cancer drivers
-   [Cancer Gene Census](https://cancer.sanger.ac.uk/census) - manually curated resource on cancer genes (soma and germline)
-   [DNA repair genes](https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html) - collection of genes involved in DNA repair
-   [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/) - collections of cancer gene panels used in clinical diagnostics

### Installation

`remotes::install_github('sigven/geneOncoXREF')`

### Usage

The package offers currently five different functions, that each retrieves a specific dataset that can be of use for gene annotation purposes.

-   `get_basic()` - retrieves basic, non-transcript-specific gene annotations. Includes tumor suppressor gene/oncogene/driver annotations from multiple resources, NCBI gene summary descriptions, as well as multiple predictions/scores when it comes to gene indispensability and loss-of-function tolerance

-   `get_gencode()` - retrieves two datasets ( *grch37* and *grch38* ) with human gene transcripts from GENCODE, including cross-references to RefSeq, [UniProt](https://www.uniprot.org), [APPRIS](https://appris.bioinfo.cnio.es/#/), and [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)

-   `get_alias()` - retrieves a list of gene synonyms, indicating which synonyms are ambiguous or nonambiguous (with respect to primary gene symbols)

-   `get_predisposition()` - retrieves a list of genes of relevance for cancer predisposition, utilizing multiple resources, including Cancer Gene Census, Genomics England PanelApp, [TCGA's PanCancer study](https://pubmed.ncbi.nlm.nih.gov/29625052/), and others.

-   `get_panels()` - retrieves a collection of > 40 different panels for various
cancer conditions, as found in Genomics England PanelApp

### Contact

sigven AT ifi.uio.no
