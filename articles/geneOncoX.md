# Getting started

  
  

## Installation

**IMPORTANT NOTE**: *geneOncoX* requires that you have R version 4.1 or
higher installed

``` r

if (!("remotes" %in% installed.packages())) {
  install.packages("remotes")
}
remotes::install_github('sigven/geneOncoX')
```

  
  

## Get basic gene annotations

This shows how to retrieve basic gene and cancer relevant gene
annotations, including how to retrieve tumor suppressor genes,
proto-oncogenes, and predicted cancer driver genes.

``` r
library(geneOncoX)

## load the data
download_dir <- tempdir()

gene_basic <- get_basic(cache_dir = download_dir)

## Number of records
nrow(gene_basic$records)
#> [1] 65334

## Show metadata for underlying resources
gene_basic$metadata
#>                             source
#> 1                       CancerMine
#> 2          Network of Cancer Genes
#> 3                          IntoGen
#> 4                             NCBI
#> 5        Bailey et al., Cell, 2018
#> 6  Sanchez-Vega et al., Cell, 2018
#> 7                            F1CDx
#> 8                           TSO500
#> 9         DNA repair gene database
#> 10                          dbNSFP
#> 11                            CPIC
#>                                                                       source_description
#> 1    Predicted tumor suppressors/oncogenes/cancer drivers from text mining of literature
#> 2                                                            Tumor suppressors/oncogenes
#> 3                                                          Predicted cancer driver genes
#> 4  Basic gene identifiers (symbol, entrez ID, HGNC ID, name, synonyms, function summary)
#> 5                     Predicted cancer driver genes / likely false positive driver genes
#> 6                                   Collection of curated signalling pathway annotations
#> 7                Collection of genes covered by Foundation One's F1CDx cancer gene panel
#> 8                     Collection of genes covered by Illumina's TSO500 cancer gene panel
#> 9                                             Collection of genes involved in DNA repair
#> 10       Gene indispensability prediction, loss-of-function intolerance predictions etc.
#> 11        Clinical Pharmacogenomics Implementation Consortium - gene/oncology drug pairs
#>                                                                   source_url
#> 1                                         http://bionlp.bcgsc.ca/cancermine/
#> 2                                           http://network-cancer-genes.org/
#> 3                                           https://www.intogen.org/download
#> 4                                         https://www.ncbi.nlm.nih.gov/gene/
#> 5                                  https://pubmed.ncbi.nlm.nih.gov/29625053/
#> 6                                  https://pubmed.ncbi.nlm.nih.gov/29625050/
#> 7                  https://www.foundationmedicine.com/test/foundationone-cdx
#> 8  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0260089
#> 9                         https://panelapp.genomicsengland.co.uk/panels/256/
#> 10                              https://sites.google.com/site/jpopgen/dbNSFP
#> 11                                          https://cpicpgx.org/genes-drugs/
#>                                            source_citation    source_version
#> 1                Lever et al., Nat Methods, 2019; 31110280 v51 (August 2025)
#> 2               Repana et al., Genome Biol, 2019; 30606230              v7.2
#> 3  Martínez-Jiménez et al., Nat Rev Cancer, 2020; 32778778        2024.09.20
#> 4          Brown et al., Nucleic Acids Res, 2015; 25355515        2026-01-02
#> 5                      Bailey et al., Cell, 2018; 29625053              <NA>
#> 6                Sanchez-Vega et al., Cell, 2018; 29625050              <NA>
#> 7                                                     <NA>              <NA>
#> 8                                                     <NA>              <NA>
#> 9                    Woods et al., Science, 2001; 11181991              v1.1
#> 10                  Liu et al., Genome Med, 2020; 33261662              v5.3
#> 11           Caudle et al.,Curr Drug Metab, 2014; 24479687              <NA>
#>    source_abbreviation                         source_license
#> 1           cancermine                                CC0 1.0
#> 2                  ncg                       Free/open access
#> 3              intogen                                CC0 1.0
#> 4                 ncbi               NCBI data usage policies
#> 5           bailey2018                       Free/open access
#> 6      sanchezvega2018                       Free/open access
#> 7       foundation_one                       Free/open access
#> 8             illumina                       Free/open access
#> 9      woods_dnarepair                       Free/open access
#> 10              dbnsfp Free for non-commercial, academic use.
#> 11                cpic                                CC0 1.0
#>                                        source_license_url
#> 1      https://creativecommons.org/publicdomain/zero/1.0/
#> 2                                                    <NA>
#> 3      https://creativecommons.org/publicdomain/zero/1.0/
#> 4  https://www.ncbi.nlm.nih.gov/home/about/policies/#data
#> 5                                                    <NA>
#> 6                                                    <NA>
#> 7                                                    <NA>
#> 8                                                    <NA>
#> 9                                                    <NA>
#> 10                                                   <NA>
#> 11                                                   <NA>
```

  
  

### Get classified tumor suppressor genes

``` r

## Get tumor suppressor genes - as indicated from either
## Network of Cancer Genes (NCG)
## - show literature support from CancerMine

tsg <- gene_basic$records |>
  dplyr::filter(ncg_tsg == TRUE) |>
  dplyr::select(symbol, entrezgene, name, gene_biotype,
                ncg_tsg, cancermine_cit_links_tsg, 
                ncbi_function_summary) |>
  dplyr::rename(support_cancermine = cancermine_cit_links_tsg)

## Make as datatable
tsg_table <- DT::datatable(
  tsg, 
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), dom = "Bfrtip"))
```

  

  
  

### Get classified proto-oncogenes

``` r

## Get proto-oncogenes - as indicated from either
## Cancer Gene Census (CGC) or Network of Cancer Genes (NCG),
## - show literature support from CancerMine
oncogene <- gene_basic$records |>
  dplyr::filter(ncg_oncogene == TRUE) |>
  dplyr::select(symbol, entrezgene, 
                name, gene_biotype,
                ncg_oncogene,
                cancermine_cit_links_oncogene, 
                ncbi_function_summary) |>
  dplyr::rename(support_cancermine = 
                  cancermine_cit_links_oncogene)

## Make as datatable
oncogene_table <- DT::datatable(
  oncogene, 
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip"))
```

  

  
  

### Get predicted cancer driver genes

``` r

## Get predicted cancer driver genes - as indicated from either
## - Network of Cancer Genes (NCG) - canonical drivers
## - IntOGen mutational driver catalogue
## - TCGA's PanCancer driver prediction (Bailey et al., Cell, 2018)
##
##  Rank hits by how many sources that contribute to classification
##
cancer_driver <- gene_basic$records |>
  dplyr::filter(intogen_driver == TRUE |
                ncg_driver == TRUE |
                tcga_driver == TRUE) |>
  dplyr::mutate(driver_score = 0) |>
  # dplyr::mutate(driver_score = dplyr::if_else(
  #   cgc_driver_tier1 == T, 
  #   driver_score + 2,
  #   as.numeric(driver_score))) |>
  # dplyr::mutate(driver_score = dplyr::if_else(
  #   cgc_driver_tier2 == T, 
  #   driver_score + 1,
  #   as.numeric(driver_score))) |>
  dplyr::mutate(driver_score = dplyr::if_else(
    intogen_driver == T, 
    driver_score + 1,
    as.numeric(driver_score))) |>
  dplyr::mutate(driver_score = dplyr::if_else(
    ncg_driver == T, 
    driver_score + 1,
    as.numeric(driver_score))) |>
  dplyr::mutate(driver_score = dplyr::if_else(
    tcga_driver == T, 
    driver_score + 1,
    as.numeric(driver_score))) |>
  dplyr::select(symbol, entrezgene, 
                name, gene_biotype,
                intogen_driver,
                intogen_role,
                ncg_driver,
                ncg_phenotype,
                tcga_driver,
                driver_score,
                cancermine_cit_links_driver, 
                ncbi_function_summary) |>
  dplyr::rename(support_cancermine = 
                  cancermine_cit_links_driver) |>
  dplyr::arrange(dplyr::desc(driver_score), symbol)

## Make as datatable
driver_table <- DT::datatable(
  cancer_driver, 
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip"))
```

  

  
  

## Get cancer predisposition genes

This show how to retrieve known cancer predisposition genes, utilizing
multiple sources, including Cancer Gene Census, Genomics England
PanelApp, TCGA’s pan-cancer study of germline variants, and
other/user-curated entries.

``` r

## load the data
gene_predisposition <- get_predisposition(cache_dir = download_dir)

## Number of cancer predisposition genes
nrow(gene_predisposition$records |> dplyr::filter(
  !stringr::str_detect(cpg_source, "^(ACMG_SF|CPIX_PGX_ONCOLOGY)$")
))
#> [1] 559

## Get statistics regarding how reference sources on 
## cancer predisposition genes contribute
##
## CGC - Cancer Gene Census (germline)
## PANEL_APP - N = 43 gene panels for inherited cancer conditions/
##             cancer syndromes (Genomics England PanelApp)
## CURATED_OTHER - curated/user-contributed genes
## TCGA_PANCAN_2018 - TCGA's pancancer analysis of
##                    germline variants in cancer
##                    (Huang et al., Cell, 2019)
plyr::count(gene_predisposition$records$cpg_source) |>
  dplyr::arrange(dplyr::desc(freq)) |>
  dplyr::filter(x != "ACMG_SF") |>
  dplyr::filter(x != "CPIC_PGX_ONCOLOGY")
#>                                               x freq
#> 1                                     PANEL_APP  275
#> 2                                 CURATED_OTHER  116
#> 3          CANVAR_UK&PANEL_APP&TCGA_PANCAN_2018   68
#> 4                    PANEL_APP&TCGA_PANCAN_2018   34
#> 5  ACMG_SF&CANVAR_UK&PANEL_APP&TCGA_PANCAN_2018   24
#> 6                           CANVAR_UK&PANEL_APP   15
#> 7                              TCGA_PANCAN_2018   12
#> 8            ACMG_SF&PANEL_APP&TCGA_PANCAN_2018    4
#> 9                             ACMG_SF&PANEL_APP    2
#> 10                                    CANVAR_UK    2
#> 11                   CANVAR_UK&TCGA_PANCAN_2018    2
#> 12           ACMG_SF&CANVAR_UK&TCGA_PANCAN_2018    1
#> 13                     ACMG_SF&TCGA_PANCAN_2018    1


## Cancer predisposition metadata
gene_predisposition$metadata
#>                                        source
#> 1                   Genomics England PanelApp
#> 2                                        NCBI
#> 3                    Huang et al., Cell, 2018
#> 4        Maxwell et al., Am J Hum Genet, 2016
#> 5 Cancer predisposition genes - curated/other
#> 6                   ACMG - secondary findings
#> 7                                   CanVar-UK
#>                                                                               source_description
#> 1 Collection of > 40 dedicated gene panels for various inherited cancer conditions and syndromes
#> 2          Basic gene identifiers (symbol, entrez ID, HGNC ID, name, synonyms, function summary)
#> 3                   Collection of cancer predisposition genes screened in TCGA's pancancer study
#> 4                                      Mechanisms of inheritance for cancer predisposition genes
#> 5                         Candidate cancer predisposition genes - contributed e.g. by CPSR users
#> 6            Genes recommended for reporting of incidental findings in clinical exome sequencing
#> 7                                                    Cancer predisposition gene variant database
#>                                  source_url
#> 1   https://panelapp.genomicsengland.co.uk/
#> 2        https://www.ncbi.nlm.nih.gov/gene/
#> 3 https://pubmed.ncbi.nlm.nih.gov/29625052/
#> 4 https://pubmed.ncbi.nlm.nih.gov/27153395/
#> 5                                      <NA>
#> 6 https://pubmed.ncbi.nlm.nih.gov/40568962/
#> 7                     https://canvaruk.org/
#>                                   source_citation source_version
#> 1        Martin et al., Nat Genet, 2019; 31676867       v1 (API)
#> 2 Brown et al., Nucleic Acids Res, 2015; 25355515     2026-01-02
#> 3              Huang et al., Cell, 2018; 29625052           <NA>
#> 4  Maxwell et al., Am J Hum Genet, 2016; 27153395           <NA>
#> 5                                            <NA>       20251106
#> 6           Lee et al., Genet Med, 2025; 40568962           v3.3
#> 7                                            <NA>           v2.2
#>   source_abbreviation                             source_license
#> 1                gepa Commercial use requires separate agreement
#> 2                ncbi                   NCBI data usage policies
#> 3    tcga_pancan_2018                           Free/open access
#> 4         maxwell2016                           Free/open access
#> 5           cpg_other                           Free/open access
#> 6             acmg_sf                           Free/open access
#> 7           canvar_uk                           Free/open access
#>                                                                                 source_license_url
#> 1 https://panelapp.genomicsengland.co.uk/media/files/GEL_-_PanelApp_Terms_of_Use_December_2019.pdf
#> 2                                           https://www.ncbi.nlm.nih.gov/home/about/policies/#data
#> 3                                                                                             <NA>
#> 4                                                                                             <NA>
#> 5                                                                                             <NA>
#> 6                                                  https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
#> 7                                                                                             <NA>
```

  
  

## Get cancer gene panels

This shows how to retrieve genes from cancer gene panels defined in
Genomics England PanelApp.

``` r

## load the data
gene_panels <- get_panels(cache_dir = download_dir)

## panel data for genome build grch38
panel_data <- gene_panels$records |>
  dplyr::filter(genome_build == "grch38")

## show number of genes in each panel
gene_freq <- as.data.frame(panel_data |>
  dplyr::group_by(gepa_panel_name) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::arrange(desc(n)) |>
  dplyr::rename(panel_name = gepa_panel_name))

gene_freq
#>                                                                  panel_name   n
#> 1                          DNA Repair Genes pertinent cancer susceptibility 178
#> 2                                                   Childhood solid tumours 122
#> 3                         Haematological malignancies cancer susceptibility 108
#> 4                                 Adult solid tumours cancer susceptibility 105
#> 5                              Haematological malignancies for rare disease  89
#> 6                             Childhood solid tumours cancer susceptibility  85
#> 7                                      Adult solid tumours for rare disease  58
#> 8                                    Multiple monogenic benign skin tumours  46
#> 9                                                    Sarcoma susceptibility  44
#> 10                                            Sarcoma cancer susceptibility  33
#> 11                                                         GI tract tumours  31
#> 12                                                 Neurofibromatosis Type 1  30
#> 13                                   Inherited non-medullary thyroid cancer  29
#> 14                                                   Familial breast cancer  27
#> 15                         Inherited ovarian cancer (without breast cancer)  27
#> 16    Familial Tumours Syndromes of the central & peripheral Nervous system  21
#> 17                            Inherited phaeochromocytoma and paraganglioma  20
#> 18 Inherited polyposis and early onset colorectal cancer - germline testing  19
#> 19                                                Familial rhabdomyosarcoma  18
#> 20                                                   Inherited renal cancer  18
#> 21                Inherited predisposition to acute myeloid leukaemia (AML)  16
#> 22                                               Multiple endocrine tumours  16
#> 23                                                 Familial prostate cancer  14
#> 24                        Colorectal cancer pertinent cancer susceptibility  13
#> 25                                         Genodermatoses with malignancies  13
#> 26                                              Inherited pancreatic cancer  13
#> 27                     Head and neck cancer pertinent cancer susceptibility  12
#> 28                    Neuroendocrine cancer pertinent cancer susceptibility  12
#> 29                             Renal cancer pertinent cancer susceptibility  10
#> 30                           Ovarian cancer pertinent cancer susceptibility   9
#> 31                                                        Familial melanoma   8
#> 32                             Brain cancer pertinent cancer susceptibility   7
#> 33                            Breast cancer pertinent cancer susceptibility   7
#> 34                                         Inherited predisposition to GIST   7
#> 35                                                       Parathyroid Cancer   7
#> 36                       Endometrial cancer pertinent cancer susceptibility   6
#> 37                                Inherited MMR deficiency (Lynch syndrome)   5
#> 38                          Prostate cancer pertinent cancer susceptibility   5
#> 39                           Thyroid cancer pertinent cancer susceptibility   5
#> 40                           Bladder cancer pertinent cancer susceptibility   4
#> 41            Upper gastrointestinal cancer pertinent cancer susceptibility   4
#> 42                                 Melanoma pertinent cancer susceptibility   3
#> 43                                                Familial rhabdoid tumours   2
#> 44         Inherited susceptibility to acute lymphoblastoid leukaemia (ALL)   2
```

  
  

## Get gene aliases

This shows how to retrieve ambiguous and unambiguous gene aliases
(i.e. with respect to primary gene symbols).

``` r

## load the data
gene_alias <- get_alias(cache_dir = download_dir)

## number of gene synonyms that are ambiguous
nrow(dplyr::filter(gene_alias$records, ambiguous == TRUE))
#> [1] 7164

## show structure of alias records
head(gene_alias$records)
#>        alias      symbol entrezgene n_primary_map ambiguous source
#> 1   'C-K-RAS        KRAS       3845             1     FALSE   NCBI
#> 2     (FM-3)       NMUR1      10316             1     FALSE   NCBI
#> 3    (IV)-44 IGHVIV-44-1      28337             1     FALSE   NCBI
#> 4      (P)RR     ATP6AP2      10159             1     FALSE   NCBI
#> 5 (ppGpp)ase       HDDC3     374659             1     FALSE   NCBI
#> 6 A-116A10.1        NAE1       8883             1     FALSE   NCBI
#>   is_primary_symbol other_index
#> 1             FALSE        <NA>
#> 2             FALSE        <NA>
#> 3             FALSE        <NA>
#> 4             FALSE        <NA>
#> 5             FALSE        <NA>
#> 6             FALSE        <NA>
```

  
  

## Get GENCODE transcripts

This shows how to retrieve GENCODE transcripts for `grch37` and
`grch38`.

``` r

## load the data
gene_gencode <- get_gencode(cache_dir = download_dir)

## number of transcript records - grch37
nrow(gene_gencode$records$grch37)
#> [1] 196520

## number of transcript records - grch38
nrow(gene_gencode$records$grch38)
#> [1] 507365

## show colnames for transcript records
colnames(gene_gencode$records$grch38)
#>  [1] "chrom"                      "start"                     
#>  [3] "end"                        "transcript_start"          
#>  [5] "transcript_end"             "cds_start"                 
#>  [7] "strand"                     "ensembl_gene_id"           
#>  [9] "ensembl_gene_id_full"       "ensembl_transcript_id"     
#> [11] "ensembl_transcript_id_full" "ensembl_protein_id"        
#> [13] "symbol"                     "symbol_gencode"            
#> [15] "hgnc_id"                    "entrezgene"                
#> [17] "name"                       "gene_biotype"              
#> [19] "transcript_biotype"         "tag"                       
#> [21] "refseq_protein_id"          "refseq_transcript_id"      
#> [23] "mane_select"                "mane_plus_clinical"        
#> [25] "refseq_select"              "principal_isoform_flag"    
#> [27] "uniprot_acc"                "uniprot_id"                
#> [29] "ensembl_version"            "gencode_version"           
#> [31] "uniprot_version"

## show metadata for underlying resources
gene_gencode$metadata
#>            source
#> 1         GENCODE
#> 2 Ensembl Biomart
#> 3       UniprotKB
#> 4          APPRIS
#>                                                         source_description
#> 1                                                   Human gene transcripts
#> 2 API for retrieval of gene and transcript cross-references (MANE, RefSeq)
#> 3      UniProt identifiers and accessions with cross-references to Ensembl
#> 4                                 Prinicipal transcript isoform annotation
#>                                       source_url
#> 1                  https://www.gencodegenes.org/
#> 2       https://www.ensembl.org/biomart/martview
#> 3                        https://www.uniprot.org
#> 4 https://apprisws.bioinfo.cnio.es/landing_page/
#>                                         source_citation source_version
#> 1    Frankish et al., Nucleic Acids Res, 2021; 33270111             49
#> 2  Cunningham et al., Nucleic Acids Res, 2022; 34791404            115
#> 3 UniProt Consortium, Nucleic Acids Res, 2021; 33237286        2025_03
#> 4    Rodriguez et al, Nucleic Acids Res, 2022; 34755885     2026-01-02
#>   source_abbreviation        source_license
#> 1             gencode      Free/open access
#> 2             ensembl EMBL-EBI terms of use
#> 3             uniprot             CC BY 4.0
#> 4              appris      Free/open access
#>                             source_license_url
#> 1                                         <NA>
#> 2     https://www.ebi.ac.uk/about/terms-of-use
#> 3 https://creativecommons.org/licenses/by/4.0/
#> 4                                         <NA>
```

  
  

## Get DNA repair genes

``` r

## load the data
gene_dna_repair <- gene_basic$records |>
  dplyr::filter(!is.na(woods_dnarepair_class)) |>
  dplyr::select(symbol, woods_dnarepair_class,
                woods_dnarepair_activity)

## count number of genes in each class
dna_repair_class_freq <- as.data.frame(
  gene_dna_repair |>
  dplyr::group_by(woods_dnarepair_class) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::arrange(desc(n)) |>
  dplyr::rename(dna_repair_class = woods_dnarepair_class))

dna_repair_class_freq
#>                                                                  dna_repair_class
#> 1                                                Nucleotide excision repair (NER)
#> 2                                                        Homologous recombination
#> 3                                                                  Fanconi anemia
#> 4                                            DNA polymerases (catalytic subunits)
#> 5                                       Other conserved DNA damage response genes
#> 6                                                      Base excision repair (BER)
#> 7                                                 Ubiquitination and modification
#> 8                                                  Mismatch excision repair (MMR)
#> 9              Other identified genes with known or suspected DNA repair function
#> 10                                               Editing and processing nucleases
#> 11                                                     Non-homologous end-joining
#> 12                                     Other BER and strand break joining factors
#> 13 Genes defective in diseases associated with sensitivity to DNA damaging agents
#> 14                                           Chromatin Structure and Modification
#> 15                                                      Direct reversal of damage
#> 16                                                 Modulation of nucleotide pools
#> 17                    Poly(ADP-ribose) polymerase (PARP) enzymes that bind to DNA
#> 18                                         Repair of DNA-topoisomerase crosslinks
#> 19                                                                        Unknown
#>     n
#> 1  28
#> 2  21
#> 3  17
#> 4  15
#> 5  15
#> 6  11
#> 7  11
#> 8  10
#> 9   9
#> 10  8
#> 11  7
#> 12  6
#> 13  5
#> 14  3
#> 15  3
#> 16  3
#> 17  3
#> 18  2
#> 19  1
```

  
  

## Get TSO500 genes

``` r

## load the data
gene_tso500 <- gene_basic$records |>
  dplyr::filter(!is.na(illumina_tso500))

## number of genes covered by TSO500 panel
nrow(gene_tso500)
#> [1] 523

## types of variant types covered
illumina_tso500_variant_freq <- as.data.frame(
  gene_tso500 |>
  dplyr::group_by(illumina_tso500) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::arrange(desc(n)))

illumina_tso500_variant_freq
#>                 illumina_tso500   n
#> 1                     SNV_INDEL 433
#> 2            CNA_GAIN,SNV_INDEL  33
#> 3          RNA_FUSION,SNV_INDEL  31
#> 4 CNA_GAIN,RNA_FUSION,SNV_INDEL  22
#> 5 CNA_LOSS,RNA_FUSION,SNV_INDEL   2
#> 6            CNA_LOSS,SNV_INDEL   2
```

  
  

## Session Info

``` r
# set eval = FALSE if you don't want this info (useful for reproducibility) 
# to appear
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] geneOncoX_1.2.7
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0    dplyr_1.1.4       compiler_4.5.2    crayon_1.5.3     
#>  [5] Rcpp_1.1.0        tidyselect_1.2.1  stringr_1.6.0     jquerylib_0.1.4  
#>  [9] systemfonts_1.3.1 textshaping_1.0.4 yaml_2.3.12       fastmap_1.2.0    
#> [13] plyr_1.8.9        R6_2.6.1          generics_0.1.4    curl_7.0.0       
#> [17] knitr_1.51        htmlwidgets_1.6.4 tibble_3.3.0      desc_1.4.3       
#> [21] bslib_0.9.0       pillar_1.11.1     rlang_1.1.6       DT_0.34.0        
#> [25] stringi_1.8.7     cachem_1.1.0      lgr_0.5.0         xfun_0.55        
#> [29] fs_1.6.6          sass_0.4.10       otel_0.2.0        cli_3.6.5        
#> [33] withr_3.0.2       pkgdown_2.2.0     magrittr_2.0.4    crosstalk_1.2.2  
#> [37] digest_0.6.39     lifecycle_1.0.4   vctrs_0.6.5       evaluate_1.0.5   
#> [41] gargle_1.6.0      glue_1.8.0        ragg_1.5.0        googledrive_2.1.2
#> [45] httr_1.4.7        rmarkdown_2.30    purrr_1.2.0       tools_4.5.2      
#> [49] pkgconfig_2.0.3   htmltools_0.5.9
```

  
  

## References

Bailey, Matthew H, Collin Tokheim, Eduard Porta-Pardo, Sohini Sengupta,
Denis Bertrand, Amila Weerasinghe, Antonio Colaprico, et al. 2018.
“Comprehensive Characterization of Cancer Driver Genes and Mutations.”
*Cell* 173 (2): 371–385.e18.
<http://dx.doi.org/10.1016/j.cell.2018.02.060>.

Lever, Jake, Eric Y Zhao, Jasleen Grewal, Martin R Jones, and Steven J M
Jones. 2019. “CancerMine: A Literature-Mined Resource for Drivers,
Oncogenes and Tumor Suppressors in Cancer.” *Nat. Methods* 16 (6):
505–7. <http://dx.doi.org/10.1038/s41592-019-0422-y>.

Liu, Xiaoming, Chang Li, Chengcheng Mou, Yibo Dong, and Yicheng Tu.
2020. “dbNSFP V4: A Comprehensive Database of Transcript-Specific
Functional Predictions and Annotations for Human Nonsynonymous and
Splice-Site SNVs.” *Genome Med.* 12 (1): 103.
<http://dx.doi.org/10.1186/s13073-020-00803-9>.

Martin, Antonio Rueda, Eleanor Williams, Rebecca E Foulger, Sarah Leigh,
Louise C Daugherty, Olivia Niblock, Ivone U S Leong, et al. 2019.
“PanelApp Crowdsources Expert Knowledge to Establish Consensus
Diagnostic Gene Panels.” *Nat. Genet.* 51 (November): 1560–65.
<http://dx.doi.org/10.1038/s41588-019-0528-2>.

Martı́nez-Jiménez, Francisco, Ferran Muiños, Inés Sentı́s, Jordi Deu-Pons,
Iker Reyes-Salazar, Claudia Arnedo-Pac, Loris Mularoni, et al. 2020. “A
Compendium of Mutational Cancer Driver Genes.” *Nat. Rev. Cancer* 20
(10): 555–72. <https://www.nature.com/articles/s41568-020-0290-x>.

Repana, Dimitra, Joel Nulsen, Lisa Dressler, Michele Bortolomeazzi,
Santhilata Kuppili Venkata, Aikaterini Tourna, Anna Yakovleva, Tommaso
Palmieri, and Francesca D Ciccarelli. 2019. “The Network of Cancer Genes
(NCG): A Comprehensive Catalogue of Known and Candidate Cancer Genes
from Cancer Sequencing Screens.” *Genome Biol.* 20 (1): 1.
<http://dx.doi.org/10.1186/s13059-018-1612-0>.

Sondka, Zbyslaw, Sally Bamford, Charlotte G Cole, Sari A Ward, Ian
Dunham, and Simon A Forbes. 2018. “The COSMIC Cancer Gene Census:
Describing Genetic Dysfunction Across All Human Cancers.” *Nat. Rev.
Cancer* 18 (11): 696–705. <http://dx.doi.org/10.1038/s41568-018-0060-1>.

Wood, R D, M Mitchell, J Sgouros, and T Lindahl. 2001. “Human DNA Repair
Genes.” *Science* 291 (5507): 1284–89.
<http://dx.doi.org/10.1126/science.1056154>.
