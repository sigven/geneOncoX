# Get human GENCODE transcript dataset

Downloads and returns a pre-processed dataset of GENCODE transcripts.
The dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a list with two data frames of transcripts, one for
  `grch37`, and one for `grch38` (see specific column format below)

The GENCODE datasets provided were established by parsing the GENCODE
[comprehensive gene annotation GTF
file](https://www.gencodegenes.org/human/) with the
[valr](https://rnabioco.github.io/valr/) package, and also providing
cross-references to RefSeq (through the use of the
[biomaRt](https://bioconductor.org/packages/biomaRt/) package), and
UniProt accessions/identifiers.

## Usage

``` r
get_gencode(
  cache_dir = NA,
  force_download = FALSE,
  chromosomes_only = TRUE,
  ensembl_release = 114
)
```

## Source

<https://www.gencodegenes.org/>

<https://uniprot.org>

<https://www.ncbi.nlm.nih.gov/gene>

<https://apprisws.bioinfo.cnio.es/landing_page/>

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

- chromosomes_only:

  Logical indicating if transcripts returned should belong to ordinary
  chromosomes only (scaffolds etc. should be ignored)

- ensembl_release:

  version of Ensembl to use - this will dictate the version of GENCODE
  used for grch38

## Value

**metadata** - A data frame with 3 rows and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used

- *abbreviation* - abbreviation (e.g. used in column names of records)

**records** - A list with one data frame per assembly

- *chrom* - chromosome

- *start* - transcript start with 5kb padding (upstream)

- *end* - transcript end with 5kb padding (downstream)

- *transcript_start* - transcript start

- *transcript_end* - transcript end

- *strand* - strand

- *ensembl_gene_id* - Ensembl gene identifier

- *ensembl_transcript_id* - Ensembl transcript identifier

- *ensembl_transcript_id_full* - Ensembl transcript identifier (with
  version)

- *ensembl_protein_id* - Ensembl protein identifier

- *symbol* - official gene symbol

- *hgnc_id* - HGNC gene identifier

- *entrezgene* - Entrez gene identifier

- *name* - genename

- *gene_biotype* - [gene
  biotype](https://www.gencodegenes.org/pages/biotypes.html) (GENCODE)

- *transcript_biotype* - [transcript
  biotype](https://www.gencodegenes.org/pages/biotypes.html) (GENCODE))

- *tag* - [tag](https://www.gencodegenes.org/pages/biotypes.html)
  (GENCODE))

- *refseq_protein_id* - RefSeq peptide identifier

- *refseq_transcript_id* - RefSeq mRNA identifier

- *mane_select* - [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)
  Select

- *mane_plus_clinical* -
  [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) Plus Clinical

- *principal_isoform_flag* - Principal [isoform
  tag](https://appris.bioinfo.cnio.es/#/help/scores) (APPRIS)

- *uniprot_acc* - UniProtKB protein accession

- *uniprot_id* - UniProtKB protein identifier

- *ensembl_version* - Ensembl version

- *gencode_version* - GENCODE version

- *uniprot_version* - UniProtKB version

list object with GENCODE transcripts (records) and associated metadata
(metadata)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geneOncoX)
transcripts_gencode <- get_gencode(cache_dir = tempdir())
} # }
```
