pcgrdb_dir <- list()
pcgrdb_dir[["data-raw"]] <- 
  "/Users/sigven/project_data/packages/package__geneOncoX/geneOncoX/data-raw/gd_local"
ensembl_release <- 112

gene_oncox <- list()
gene_oncox[['basic']] <-
  geneOncoX::get_basic(
    cache_dir = pcgrdb_dir[['data-raw']])

##  Get TSG/oncogene/cancer driver support per gene -
##  utilizing multiple sources
##   - Cancer Gene Census (CGC)
##   - Network of Cancer Genes (NCG)
##   - CancerMine
##   - IntOGen
##   - TCGA
cancer_gene_support <-
  geneOncoX:::assign_cancer_gene_roles(
    gox_basic = gene_oncox[['basic']]) |>
  dplyr::left_join(
    dplyr::select(
      gene_oncox$basic$records,
      entrezgene, gene_biotype))

gene_oncox[['gencode']] <-
  geneOncoX::get_gencode(
    cache_dir = pcgrdb_dir[['data-raw']],
    ensembl_release = ensembl_release,
    chromosomes_only = TRUE)

gene_oncox[['gencode_all']] <-
  geneOncoX::get_gencode(
    cache_dir = pcgrdb_dir[['data-raw']],
    ensembl_release = ensembl_release,
    chromosomes_only = FALSE)


gene_oncox[['alias']] <-
  geneOncoX::get_alias(
    cache_dir = pcgrdb_dir[['data-raw']])

gene_oncox[['panels']] <-
  geneOncoX::get_panels(
    cache_dir = pcgrdb_dir[['data-raw']])

gene_oncox[['otp_rank']] <-
  geneOncoX::get_otp_rank(
    cache_dir = pcgrdb_dir[['data-raw']]
  )


biomarkers <-
  pharmOncoX:::get_biomarkers(
    cache_dir = pcgrdb_dir[['data-raw']]
  )

actionable_genes <-
  data.frame(
    'entrezgene' = unique(c(
      biomarkers$data$civic$variant$entrezgene,
      biomarkers$data$cgi$variant$entrezgene))
  ) |>
  dplyr::filter(!is.na(entrezgene)) |>
  dplyr::mutate(actionable_gene = T)

gene_oncox[['predisposition']] <-
  geneOncoX::get_predisposition(
    cache_dir = pcgrdb_dir[['data-raw']])


gene_oncox[['predisposition']]$records <-
  gene_oncox[['predisposition']]$records |>
  dplyr::mutate(
    cpg_phenotypes = dplyr::if_else(
      (nchar(cpg_phenotypes) == 0 |
         is.na(cpg_phenotypes)),
      ".",
      as.character(cpg_phenotypes)
    )
  )

gene_oncox[['panels']]$records <-
  gene_oncox[['panels']]$records |>
  dplyr::left_join(
    dplyr::select(
      gene_oncox[['predisposition']]$records,
      entrezgene, cpg_mod
    ), by = "entrezgene"
  ) |>
  dplyr::rename(
    gepa_mod = cpg_mod
  )


gene_metadata_all <- gene_oncox[['predisposition']][['metadata']] |>
  dplyr::bind_rows(gene_oncox[['panels']][['metadata']]) |>
  dplyr::bind_rows(gene_oncox[['basic']][['metadata']]) |>
  dplyr::bind_rows(gene_oncox[['gencode']][['metadata']]) |>
  dplyr::distinct()

genome_builds <- c("grch37","grch38")

for (build in genome_builds) {
  
  gene_cpg <- gene_oncox[['predisposition']]$records |>
    dplyr::left_join(
      dplyr::select(
        gene_oncox[['gencode']][['records']][[build]],
        c("entrezgene","ensembl_gene_id","gene_biotype")),
      by = c("entrezgene","gene_biotype")
    ) |>
    dplyr::distinct()
  
  gencode <-
    gene_oncox[['gencode']][['records']][[build]] |>
    dplyr::mutate(chrom = stringr::str_replace(
      chrom, "^chr", ""
    )) |>
    dplyr::left_join(
      dplyr::select(
        gene_oncox[['basic']][['records']],
        -c(hgnc_id, name,
           symbol)),
      by = c("entrezgene","gene_biotype")) |>
    dplyr::mutate(cgc_tier = dplyr::if_else(
      cgc_tier == "NA",
      as.numeric(NA),
      as.numeric(cgc_tier)
    )) |>
    dplyr::left_join(
      cancer_gene_support,
      by = c("entrezgene","gene_biotype")
    )
  
  num_tsg <- gencode |>
    dplyr::filter(tsg == TRUE) |>
    dplyr::select(entrezgene) |>
    dplyr::distinct() |>
    nrow()
  
  num_oncogene <- gencode |>
    dplyr::filter(oncogene == TRUE) |>
    dplyr::select(entrezgene) |>
    dplyr::distinct() |>
    nrow()
  
  num_cpg <- gene_cpg |>
    dplyr::select(entrezgene) |>
    dplyr::distinct() |>
    nrow()
  
  cat("Number of TSGs in", build, ":", num_tsg, "\n")
  cat("Number of oncogenes in", build, ":", num_oncogene, "\n")
  cat("Number of predisposition genes in", build, ":", num_cpg, "\n")
}