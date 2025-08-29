# Version 1.1.8

- Bug in mechanism-of-disease fixed

# Version 1.1.7

- Fix missing MANE Select and MANE Plus Clinical annotations

# Version 1.1.6

- Updated CancerMine (v51, August 2025)
- Curated information on mechanism of disease for LZTR1, RTEL1, CTC1

# Version 1.1.5

- Updated NCBI gene info (2024-08-24)
- Added information on mechanism of inheritance for more cancer 
predisposition genes (from ClinGen)

# Version 1.1.2

- GENCODE v48 (to be used for VEP v114)
- Cancer relevance ranks now based on data from OTP v25.03
- Removed some user-contributed cancer predisposition genes - no/very weak evidence for cancer relevance
  - SERPINA1, GJB2, DHCR7, ASPM, UROD, SLC25A13
- HFE now only listed in ACMG secondary findings list (not CPG candidate)
- Updated NCBI gene info

# Version 1.1.0

- Added SAMHD1 and FAM111A as susceptibility genes
- Updated NCBI gene

# Version 1.0.12

- Added MYO3A, MYO5B, CHD1L as susceptibility genes
- Updated NCBI gene, UniProt KB 2025-01

# Version 1.0.10

- Added PAH as putative CPG gene

# Version 1.0.9

- added DPYD, NUDT15 and TPMT as PgX-related genes with respect to chemotherapy 
(included among CPG genes)

# Version 1.0.7

- added custom gene aliases (KRAS, BRAF)
- Updated NCBI summary information

# Version 1.0.6

- Updated NCBI summary information
- added custom gene aliases

# Version 1.0.4

- Update curated cancer predisposition genes

# Version 1.0.2

- Cancer Gene Census v101

# Version 1.0.0

- Updated IntOGen (v2024.10)
- Updated NCBI summary information
- Updated GENCODE to v47
- moved from bump2version (deprecated) to 
[bump-my-version](https://github.com/callowayproject/bump-my-version)

# Version 0.9.9

- Updated NCBI gene information
- Added `cds_start` to GENCODE transcripts

# Version 0.9.7

- fixed bug in Entrez gene identifiers (grch37)
- fixed bug in HGNC id map (grch38)
- added missing Entrez gene identifiers (grch37/grch38)
- added missing names for genes (grch37/grch38)

# Version 0.9.2

- GENCODE transcript tracks:
   - Pull in more up-to-date MANE select transcripts from NCBI
   - Pull in RefSeq Select identifiers (grch37)
   - Use exact version numbers for transcripts/MANE select.

# Version 0.9.1

- Use up-to-date gene symbols from NCBI in GENCODE tracks, ignore outdated Entrez
gene identifiers from BioMart

# Version 0.9.0

- Update NCBI gene information

# Version 0.8.10

- Support for GENCODE v46

# Version 0.8.9

- CGC v100
- make gene biotype more uniform across GENCODE builds

# Version 0.8.6

- Added more gene aliases (one edit distance away from existing aliases, found 
as substrings in other genename designations)

# Version 0.8.5

- Including references to non-coding RNA identifiers from RefSeq in
`refseq_transcript_id` column

# Version 0.8.3

- Updated cancer gene ranks from OTP (v24.03)
- Added GENCODE entries for transcripts originating on
scaffolds etc. (not just chromosomes)

# Version 0.8.1

- Updated NCBI gene info

# Version 0.8.0

- Added cancer gene ranks from OTP

# Version 0.7.14

- Added support for GENCODE v45

# Version 0.7.13

- Updated gene descriptions, and updated helper function `assign_cancer_gene_roles`

# Version 0.7.12

- New panel from PanelApp included in cancer panels and included in list of cancer predisposition genes: _Inherited susceptibility to acute lymphoblastoid leukaemia (ALL) (GEP)_
- Updated Cancer Gene Census to v99

# Version 0.7.11

* Renamed columns of GENCODE data frame: 
     - `refseq_mrna` -> `refseq_transcript_id` 
     - `refseq_peptide` -> `refseq_protein_id`

# Version 0.7.10

* Removed trailing lines from Genomics England PanelApp
phenotype descriptions

# Version 0.7.9

* NCBI gene function update
* Helper function `get_cancer_gene_roles` updated

# Version 0.7.8

* Grab gene names/descriptions from Ensembl if missing from NCBI

# Version 0.7.6

* Allowing the possibility to query and retrieve different versions of GENCODE transcripts (using Ensembl version 110 to 105, effectively querying GENCODE v44 to GENCODE v39). This will only
effect grch38, for grch37, GENCODE v19 will always be used.

# Version 0.7.5

* Updated metadata table

# Version 0.7.4

* Updated IntOGen
* Fixed hyphens in URLs (CGC, NCG, IntOGen)
* Added licensing terms

# Version 0.7.3

* Minor change to column names in cancer predisposition gene records

# Version 0.7.2

* Updated Cancer Gene Census (v98)
* Updated gene annotations from dbNSFP (v4.4)
* Updated Network of Cancer Genes (v7.1)
* Updated UniProt (2023_02)

# Version 0.7.1

* Retrieved new PanelApp genes (March 2023)

# Version 0.7.0

* Upgraded CancerMine (v50, March 2023)

# Version 0.6.9

* Implemented new helper function `assign_cancer_gene_roles`
   - adds confidence support for tumor supressive/oncogenic roles of genes
     based on amount of literature evidence and manual curation status

# Version 0.6.8

* Upgraded GENCODE (v43), Ensembl v109

# Version 0.6.7

* Fixed a bug in the helper function `assign_cancer_gene_evidence`

# Version 0.6.6

* Upgraded CancerMine (v49, January 2023)
* Added a few custom gene aliases

# Version 0.6.5

* Ensured both primary symbols and aliases were present
in the gene alias data frame

# Version 0.6.4

* Use tidyverse code style
* Ignore genes from NCBI with biotype `biological-region`

# Version 0.6.3
  
* Improved UniProt identifier mappings for Ensembl transcripts
associated with the _grch37_ genome assembly
* Updated CGC (v97)

# Version 0.6.2

* Added a few curated cancer predisposition genes (contributed by UMCCR)
   * [CSDE1](https://www.ncbi.nlm.nih.gov/gene/7812), 
   [EGLN1](https://www.ncbi.nlm.nih.gov/gene/54583), 
   and [EGLN2](https://www.ncbi.nlm.nih.gov/gene/112398)
* Updates from NCBI gene info (2022-11-28)
* Updated cross-references from UniProt KB (2022_04)

# Version 0.6.1

* Updated GENCODE (v42), and NCBI gene info (2022-10-31)

# Version 0.6.0

* Added a `NEWS.md` file to track changes to the package.
* Added `geneOncoX.Rmd` to vignettes - getting started guide
* Updated code and documentation according to style
recommendations

# Version 0.5.5

* Added F1CDx annotation (Foundation One)

