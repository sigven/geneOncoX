get_panel_app_genes <-
  function(gene_info = NULL,
           build = "grch37") {
    lgr::lgr$info(
      paste0(
        "Retrieving Genomics England virtual gene panels - ",
        "cancer phenotypes"
      )
    )

    panel_id_file <-
      file.path(
        "data-raw",
        "predisposition",
        "ge_panelapp",
        "all_cancer_panels.txt"
      )
    panels <- read.table(
      file = panel_id_file, header = TRUE, sep = "\t",
      stringsAsFactors = FALSE, quote = ""
    )
    i <- 1
    all_panels <- data.frame()
    while (i <= nrow(panels)) {
      api_url <-
        paste0(
          "https://panelapp.genomicsengland.co.uk/api/v1/panels/",
          panels[i, ]$id
        )
      panel_name <- panels[i, ]$name
      d <- jsonlite::fromJSON(api_url)

      gene_phenotypes <- c()
      for (n in 1:length(d$genes$phenotypes)) {
        phenotypes <- paste(
          unique(stringr::str_replace_all(
            d$genes$phenotypes[[n]], "\\{|\\}|\\?", ""
          )),
          collapse = ";"
        )
        gene_phenotypes <- c(gene_phenotypes, phenotypes)
      }

      if (build == "grch38") {
        df <- data.frame(
          "hgnc_id" =
            stringr::str_replace(
              d$genes$gene_data$hgnc_id,
              "HGNC:", ""
            ),
          "ensembl_gene_id" =
            d$genes$gene_data$ensembl_genes$GRch38$`90`$ensembl_id,
          "gepa_moi" = d$genes$mode_of_inheritance,
          "gepa_penetrance" = d$genes$penetrance,
          "gepa_confidence_level" =
            as.integer(d$genes$confidence_level),
          "gepa_panel_name" = d$name,
          "gepa_panel_id" = d$id,
          "gepa_panel_version" = d$version,
          "gepa_phenotype" = gene_phenotypes,
          stringsAsFactors = FALSE
        )
      } else {
        df <- data.frame(
          "hgnc_id" =
            stringr::str_replace(
              d$genes$gene_data$hgnc_id,
              "HGNC:", ""
            ),
          "ensembl_gene_id" =
            d$genes$gene_data$ensembl_genes$GRch37$`82`$ensembl_id,
          "gepa_moi" = d$genes$mode_of_inheritance,
          "gepa_penetrance" = d$genes$penetrance,
          "gepa_confidence_level" =
            as.integer(d$genes$confidence_level),
          "gepa_panel_name" = d$name,
          "gepa_panel_id" = d$id,
          "gepa_panel_version" = d$version,
          "gepa_phenotype" = gene_phenotypes,
          stringsAsFactors = FALSE
        )
      }
      df$id <- i
      df$gepa_panel_url <-
        paste0(
          "https://panelapp.genomicsengland.co.uk/panels/",
          df$gepa_panel_id
        )
      all_panels <- dplyr::bind_rows(all_panels, df)
      i <- i + 1
    }


    all_repair_genes <- as.data.frame(readr::read_tsv(
      file.path(
        "data-raw",
        "predisposition",
        "ge_panelapp",
        "dna_repair.tsv"
      ),
      show_col_types = F
    ) |>
      janitor::clean_names())

    repair_genes_panel <- data.frame()
    if (build == "grch38") {
      repair_genes_panel <- all_repair_genes |>
        dplyr::rename(gepa_panel_name = level4) |>
        dplyr::mutate(
          id = i, gepa_moi = NA, gepa_panel_id = 256,
          gepa_panel_version = as.character(version),
          gepa_penetrance = NA,
          gepa_phenotype = NA,
          gepa_panel_url =
            "https://panelapp.genomicsengland.co.uk/panels/256/",
          gepa_confidence_level = 0,
          hgnc_id = stringr::str_replace(hgnc, "HGNC:", ""),
          ensembl_gene_id = ensembl_id_g_rch38
        ) |>
        dplyr::select(
          id, ensembl_gene_id, hgnc_id,
          gepa_panel_name, gepa_moi, gepa_panel_id,
          gepa_penetrance, gepa_panel_version,
          gepa_panel_url, gepa_confidence_level
        )
    } else {
      repair_genes_panel <- all_repair_genes |>
        dplyr::rename(gepa_panel_name = level4) |>
        dplyr::mutate(
          id = i, 
          gepa_moi = NA, 
          gepa_panel_id = 256,
          gepa_panel_version = as.character(version),
          gepa_penetrance = NA,
          gepa_phenotype = NA,
          gepa_panel_url =
            "https://panelapp.genomicsengland.co.uk/panels/256/",
          gepa_confidence_level = 0,
          hgnc_id = stringr::str_replace(hgnc, "HGNC:", ""),
          ensembl_gene_id = ensembl_id_g_rch37
        ) |>
        dplyr::select(
          id, ensembl_gene_id, hgnc_id, gepa_phenotype,
          gepa_panel_name, gepa_moi, gepa_panel_id,
          gepa_penetrance, gepa_panel_version,
          gepa_panel_url, gepa_confidence_level
        )
    }

    all_panels <- all_panels |>
      dplyr::filter(ensembl_gene_id != "ENSG00000277027") |>
      dplyr::filter(ensembl_gene_id != "ENSG00000277925") |>
      dplyr::filter(ensembl_gene_id != "ENSG00000278815") |>
      dplyr::mutate(
        gepa_moi2 = dplyr::if_else(
          stringr::str_detect(
            gepa_moi,
            paste0(
              "MONOALLELIC, autosomal or ",
              "pseudoautosomal, NOT imprinted"
            )
          ) |
            stringr::str_detect(
              gepa_moi,
              paste0(
                "MONOALLELIC, autosomal or ",
                "pseudoautosomal, imprinted status unknown"
              )
            ),
          "AD", as.character(NA)
        )
      ) |>
      dplyr::mutate(gepa_moi2 = dplyr::if_else(
        stringr::str_detect(
          gepa_moi, "BIALLELIC, autosomal or pseudoautosomal"
        ),
        "AR", as.character(gepa_moi2)
      )) |>
      dplyr::mutate(gepa_moi2 = dplyr::if_else(
        stringr::str_detect(
          gepa_moi, "BOTH monoallelic and biallelic"
        ),
        "AD/AR", as.character(gepa_moi2)
      )) |>
      dplyr::bind_rows(repair_genes_panel) |>
      dplyr::mutate(genome_build = build) |>
      dplyr::select(-gepa_moi) |>
      dplyr::rename(gepa_moi = gepa_moi2) |>
      # dplyr::select(-gepa_moi) |>
      dplyr::left_join(
        dplyr::select(
          gene_info, hgnc_id,
          entrezgene, name
        ),
        by = "hgnc_id", multiple = "all"
      ) |>
      dplyr::rename(genename = name) |>
      dplyr::select(-hgnc_id) |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::select(
        genome_build, id, entrezgene, genename, ensembl_gene_id,
        dplyr::everything()
      ) |>
      dplyr::mutate(
        gepa_phenotype = stringr::str_replace_all(
          gepa_phenotype, 
          "[\r\n\t]", 
          "")
      )


    return(all_panels)
  }

get_acmg_secondary_findings <- function(gene_info = NULL,
                                        cache_dir = NA) {
  ontology_maps <- phenOncoX::get_aux_maps(
    cache_dir = cache_dir
  )

  umls_concept <- ontology_maps$records$umls$concept

  acmg_sf_list <- as.data.frame(
    openxlsx::read.xlsx(
      file.path(
        "data-raw", "predisposition",
        "acmg_secondary_findings",
        "acmg_secondary_findings_v3.2.xlsx"
      ),
      sheet = 1, startRow = 3
    ) |>
      janitor::clean_names() |>
      dplyr::rename(
        symbol = gene, disease_phenotype = disease_phentyope
      ) |>
      dplyr::filter(!is.na(disease_phenotype)) |>
      dplyr::mutate(variants_to_report = stringr::str_replace_all(
        variants_to_report, "All |and |\\(.+\\)", ""
      )) |>
      dplyr::mutate(disease_phenotype = stringr::str_replace_all(
        disease_phenotype, " \\(.+\\)", ""
      )) |>
      dplyr::mutate(inheritance = stringr::str_trim(inheritance)) |>
      dplyr::mutate(
        variants_to_report = stringr::str_trim(variants_to_report)
      ) |>
      dplyr::mutate(disease_phenotype = dplyr::if_else(
        stringr::str_detect(disease_phenotype, "Long-QT"),
        stringr::str_replace_all(
          disease_phenotype, " type |-", " "
        ),
        as.character(disease_phenotype)
      )) |>
      dplyr::mutate(disease_phenotype = dplyr::if_else(
        stringr::str_detect(
          disease_phenotype, "Familial medullary"
        ),
        stringr::str_replace_all(
          disease_phenotype, " cancer", " carcinoma"
        ),
        as.character(disease_phenotype)
      )) |>
      dplyr::mutate(disease_phenotype = dplyr::if_else(
        stringr::str_detect(disease_phenotype, "WT1-related "),
        as.character("Wilms tumor 1"),
        as.character(disease_phenotype)
      )) |>
      dplyr::left_join(umls_concept,
        by = c("disease_phenotype" = "cui_name"), 
        multiple = "all", relationship = "many-to-many"
      ) |>
      dplyr::select(-c(source, main_term)) |>
      dplyr::mutate(cui = dplyr::if_else(
        disease_phenotype ==
          "Hereditary paraganglioma-pheochromocytoma syndrome",
        "C1708353",
        as.character(cui)
      )) |>
      dplyr::distinct() |>
      dplyr::group_by_at(dplyr::vars(-cui)) |>
      dplyr::summarise(
        cui = paste(cui, collapse = ", "),
        .groups = "drop"
      ) |>
      # dplyr::rename(predisp_syndrome_cui = cui) |>
      dplyr::left_join(
        dplyr::select(
          gene_info, gene_biotype, symbol, entrezgene),
        by = c("symbol"), multiple = "all",
        relationship = "many-to-many"
      ) |>
      dplyr::group_by(
        entrezgene, gene_biotype, gene_mim,
        variants_to_report, inheritance, sf_list_version
      ) |>
      dplyr::summarise(
        cui = paste(sort(unique(cui)), collapse = ", "),
        disease_phenotype = paste(
          sort(unique(disease_phenotype)),
          collapse = ";"
        ),
        disorder_mim = paste(
          sort(unique(disorder_mim)),
          collapse = ";"
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(source = "ACMG_SF")
  )

  return(acmg_sf_list)
}

get_predisposition_genes_huang018 <- function(gene_info = NULL) {
  tcga_pancan2018_genes <-
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "predisposition",
        "huang2018_cell",
        "huang_pancancer_cell_2018.xlsx"
      ),
      sheet = 1
    ) |>
    janitor::clean_names() |>
    dplyr::select(
      gene_symbol, mode_of_inheritance,
      major_associated_tumor_types, cancer_syndrome_s
    ) |>
    dplyr::rename(
      symbol = gene_symbol,
      tumor_type = major_associated_tumor_types,
      cancer_syndrome = cancer_syndrome_s
    ) |>
    dplyr::mutate(
      mode_of_inheritance = stringr::str_replace(
        stringr::str_trim(mode_of_inheritance), "\n", ""
      ),
      symbol = stringr::str_trim(symbol)
    ) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      stringr::str_trim(tumor_type),
      "(Biallelic|Monoallelic) mutations:", ""
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace_all(
      stringr::str_to_title(tumor_type), "\\s?\\r\\n", ", "
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "^\\s{0,},\\s{0,}", ""
    )) |>
    dplyr::mutate(tumor_type = stringr::str_trim(tumor_type)) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "Mds", "MDS"
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "(-| )All", " ALL"
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "Bcc", "BCC"
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "Mpn", "MPN"
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "Jmml", "JMML"
    )) |>
    dplyr::mutate(tumor_type = stringr::str_replace(
      tumor_type, "Aml", "AML"
    )) |>
    dplyr::mutate(moi = dplyr::if_else(
      stringr::str_detect(
        mode_of_inheritance, "autosomal dominant"
      ),
      "AD", as.character(NA)
    )) |>
    dplyr::mutate(moi = dplyr::if_else(
      stringr::str_detect
      (mode_of_inheritance, "autosomal recessive"),
      "AR", as.character(moi)
    )) |>
    dplyr::mutate(moi = dplyr::if_else(
      stringr::str_detect(
        mode_of_inheritance, "autosomal dominant"
      ) &
        stringr::str_detect(
          mode_of_inheritance,
          "autosomal recessive"
        ), "AD/AR",
      as.character(moi)
    )) |>
    dplyr::mutate(phenotypes = paste(
      tumor_type, cancer_syndrome,
      sep = "; "
    )) |>
    dplyr::select(
      -c(cancer_syndrome, mode_of_inheritance, tumor_type)
    ) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, 
                    entrezgene, gene_biotype),
      by = c("symbol"), 
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(source = "TCGA_PANCAN_2018") |>
    # dplyr::filter(
    #   !stringr::str_detect(
    #     symbol,"^(HFE|SLC25A13)$")
    # ) |>
    dplyr::select(-symbol)
  # dplyr::mutate(moi = NA)

  tcga_2018_predisposition_umls <-
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "predisposition",
        "huang2018_cell",
        "huang_pancancer_cell_2018.xlsx"
      ),
      sheet = 4
    ) |>
    dplyr::mutate(
      predisp_syndrome_cui = stringr::str_trim(cui_syndrome),
      susceptibility_cui = stringr::str_trim(cui_susceptibility)
    ) |>
    dplyr::mutate(susceptibility_cui = stringr::str_replace_all(
      susceptibility_cui, ";", "&"
    )) |>
    dplyr::mutate(predisp_syndrome_cui = stringr::str_replace_all(
      predisp_syndrome_cui, ";", "&"
    )) |>
    dplyr::select(symbol, predisp_syndrome_cui, susceptibility_cui) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene),
      by = c("symbol"), multiple = "all"
    ) |>
    dplyr::select(-symbol) |>
    dplyr::distinct()

  tcga_pancan2018_genes <- tcga_pancan2018_genes |>
    dplyr::left_join(tcga_2018_predisposition_umls, 
                     by = "entrezgene", multiple = "all") |>
    dplyr::select(entrezgene, source, dplyr::everything())

  return(tcga_pancan2018_genes)
}

get_canvar_cpgs <- function(gene_info = NULL){
  
  cpgs <-
    readr::read_tsv(
      file.path(
        "data-raw",
        "predisposition",
        "canvar_uk",
        "cancer_predisposition_genes_20250521.txt"
      ),
      show_col_types = F
    ) |>
    janitor::clean_names() |>
    dplyr::mutate(source = "CANVAR_UK") |>
    dplyr::select(symbol, source) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        entrezgene,
        gene_biotype,
        symbol
      ),
      by = c("symbol"), multiple = "all"
    ) |>
    dplyr::select(
      entrezgene, 
      gene_biotype, 
      source
    )
  
}

get_curated_predisposition_genes <- function(gene_info = NULL) {
  curated_other <-
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "predisposition",
        "curated_genes",
        "predisposition_curated_other.xlsx"
      ),
      sheet = 1
    ) |>
    janitor::clean_names() |>
    dplyr::mutate(
      ensembl_gene_id = stringr::str_trim(ensembl_gene_id),
      source = "CURATED_OTHER"
    ) |>
    dplyr::select(ensembl_gene_id, moi, source, reference) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        entrezgene,
        gene_biotype,
        ensembl_gene_id
      ),
      by = c("ensembl_gene_id"), multiple = "all"
    ) |>
    dplyr::select(
      entrezgene, 
      gene_biotype, 
      source,
      moi, reference
    )


  curated_ncgc <-
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "predisposition",
        "curated_genes",
        "predisposition_curated_other.xlsx"
      ),
      sheet = 2
    ) |>
    janitor::clean_names() |>
    dplyr::mutate(
      symbol = stringr::str_trim(symbol),
      source = "CURATED_OTHER"
    ) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        gene_biotype,
        symbol, 
        entrezgene
      ),
      by = c("symbol"), multiple = "all"
    ) |>
    dplyr::select(
      entrezgene, 
      source, 
      gene_biotype)

  curated_predisposition_genes <- as.data.frame(
    curated_other |>
      dplyr::bind_rows(curated_ncgc) |>
      dplyr::group_by(entrezgene, source, gene_biotype) |>
      dplyr::summarise(
        moi = paste(unique(moi), collapse = ","),
        reference = paste(unique(reference), collapse = ","),
        .groups = "drop"
      ) |>
      dplyr::mutate(moi = dplyr::if_else(
        moi == "NA", as.character(NA), as.character(moi)
      )) |>
      dplyr::mutate(reference = dplyr::if_else(
        moi == "NA", as.character(NA), as.character(reference)
      ))
  )

  return(curated_predisposition_genes)
}

get_moi_clingen <- function(gene_info = NULL,
                            cancer_phenotypes = FALSE){
  
  
  clingen_gene_curation <- 
    readr::read_tsv(
      file = file.path(
        "data-raw",
        "predisposition",
        "clingen",
        "clingen_gene_curation_list.tsv"
      ),
      show_col_types = F, skip = 5, guess_max = 3000
    ) |>
    janitor::clean_names() |>
    dplyr::select(
      gene_id,
      haploinsufficiency_description,
      haploinsufficiency_score,
    ) |>
    dplyr::rename(entrezgene = gene_id)
  
  clingene_disease_validity <- 
    readr::read_csv(
      file = file.path(
        "data-raw",
        "predisposition",
        "clingen",
        "clingen_gene_disease_summary.csv"), 
        show_col_types = F, skip = 4) |> 
        janitor::clean_names() |> 
        dplyr::filter(!startsWith(gene_symbol, "++")) |> 
        dplyr::rename(moi_clingen = moi, symbol = gene_symbol) |>
    
    dplyr::rename(hgnc_id = gene_id_hgnc) |> 
    dplyr::group_by(symbol, hgnc_id) |> 
    dplyr::summarise(
      disease_label = paste(unique(disease_label), collapse=";"), 
      disease_id_mondo = paste(unique(disease_id_mondo), collapse=";"), 
      moi_clingen = paste(
        moi_clingen = paste(unique(moi_clingen), collapse=";")), 
      classification = paste(unique(sort(classification)), collapse=";"),
      .groups = "drop") |>
    dplyr::ungroup() |>
    dplyr::mutate(moi_clingen = dplyr::if_else(
      moi_clingen == "AD;AR" | moi_clingen == "AR;AD",
      "AD/AR", as.character(moi_clingen))
    ) |>
    dplyr::filter(
      classification != "Refuted" &
        classification != "No Known Disease Relationship" &
        classification != "No Known Disease Relationship;Refuted") |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        symbol, entrezgene
      ),
      by = c("symbol"), relationship = "many-to-many"
    ) |>
    dplyr::filter(!is.na(entrezgene)) |>
    dplyr::select(symbol, entrezgene, disease_label, 
                  disease_id_mondo, moi_clingen) |>
    dplyr::distinct () |>
    dplyr::left_join(clingen_gene_curation,
      by = "entrezgene", 
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      mod_clingen = dplyr::if_else(
        haploinsufficiency_score > 1,
        "LoF", 
        as.character(NA)
      )
    ) |>
    dplyr::select(-c("symbol")) |>
    dplyr::distinct()
  
  if(cancer_phenotypes == T){
    clingene_disease_validity <- clingene_disease_validity |>
      dplyr::filter(
        stringr::str_detect(
          tolower(disease_label), 
          paste0(
            "tumor|cancer|neoplasm|carcinoma|sarcoma|leukemia|lymphoma|",
            "melanoma|medulloblastoma|neuroblastoma|myeloma|glioma|",
            "cholangiocar|lynch|fraumeni|dicer|neurofibromatosis|",
            "retinoblastoma|schwannoma")
          )
        )
  }
  
  return(clingene_disease_validity)
  
}

get_moi_mod_maxwell2016 <- function(gene_info = NULL) {
  mod_moi_predisposition <-
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "predisposition",
        "maxwell2016_ajhg",
        "maxwell_ajhg_2016.xlsx"
      ),
      sheet = 1, startRow = 2
    ) |>
    janitor::clean_names() |>
    dplyr::filter(!stringr::str_detect(
      gene_name, "^Autosomal|Misc|Genes"
    )) |>
    dplyr::mutate(
      symbol = gene_name,
      mod_maxwell = lof_known_mech_of_disease_for_pvs1,
      moi_maxwell = stringr::str_trim(mode_of_inheritance_for_pm3_bp2)
    ) |>
    dplyr::mutate(mod_maxwell = dplyr::if_else(
      mod_maxwell == "YES" & !is.na(mod_maxwell),
      "LoF",
      as.character(mod_maxwell)
    )) |>
    dplyr::mutate(mod_maxwell = dplyr::if_else(
      !is.na(mod_maxwell) &
        stringr::str_detect(mod_maxwell, "NO-GOF") &
        stringr::str_detect(mod_maxwell, "YES"),
      "GoF/LoF",
      as.character(mod_maxwell)
    )) |>
    dplyr::mutate(mod_maxwell = dplyr::if_else(
      mod_maxwell == "NO-GOF" & !is.na(mod_maxwell),
      "GoF",
      as.character(mod_maxwell)
    )) |>
    dplyr::mutate(mod_maxwell = dplyr::if_else(
      mod_maxwell == "n/a" & !is.na(mod_maxwell),
      as.character(NA),
      as.character(mod_maxwell)
    )) |>
    dplyr::mutate(moi_maxwell = dplyr::if_else(
      stringr::str_detect(moi_maxwell, "n/a|AD/XLR|Mosaic|XLR"),
      as.character(NA),
      as.character(moi_maxwell)
    )) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        symbol, entrezgene
      ),
      by = c("symbol"), relationship = "many-to-many"
    ) |>
    dplyr::select(entrezgene, mod_maxwell, moi_maxwell) |>
    dplyr::distinct()

  return(mod_moi_predisposition)
}

get_predisposition_genes <- function(gene_info = NULL,
                                     build = NULL,
                                     gene_panels = NULL,
                                     cache_dir = NA) {
  cpg_collections <- list()
  
  cpg_collections[["TCGA_PANCAN_2018"]] <-
    get_predisposition_genes_huang018(gene_info = gene_info)
  
  cpg_collections[['CANVAR_UK']] <-
    get_canvar_cpgs(gene_info = gene_info)
  
  cpg_collections[["ACMG_SF"]] <-
    get_acmg_secondary_findings(
      gene_info = gene_info,
      cache_dir = cache_dir
    ) |>
    
    dplyr::select(entrezgene, gene_biotype, 
                  inheritance, disease_phenotype) |>
    dplyr::rename(moi = inheritance, 
                  phenotypes = disease_phenotype) |>
    dplyr::mutate(source = "ACMG_SF")
  cpg_collections[["CPIC_PGX_ONCOLOGY"]] <-
    get_cpic_genes(gene_info = gene_info) |>
    dplyr::filter(entrezgene == 1806 |
                     entrezgene ==  7172 |
                     entrezgene == 55270) |> ## DPYD, TPMT and NUDT15 for now
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, gene_biotype),
      by = "entrezgene") |>
    dplyr::select(-cpic_pgx_oncology) |>
    dplyr::mutate(source = "CPIC_PGX_ONCOLOGY")
  cpg_collections[["CURATED_OTHER"]] <-
    get_curated_predisposition_genes(gene_info = gene_info) |>
    dplyr::select(-reference)
  # cpg_collections[["CGC"]] <- get_cancer_gene_census(origin = "germline") |>
  #   dplyr::rename(moi = cgc_moi, 
  #                 phenotypes = cgc_phenotype_germline) |>
  #   dplyr::select(entrezgene, 
  #                 moi, phenotypes) |>
  #   dplyr::mutate(source = "CGC") |>
  #   dplyr::left_join(
  #     dplyr::select(
  #       gene_info, entrezgene, 
  #       gene_biotype
  #     ),
  #     by = "entrezgene"
  #   )

  mod_moi_predisposition <- get_moi_mod_maxwell2016(
    gene_info = gene_info)
  
  clingen_moi <- get_moi_clingen(
    gene_info = gene_info
  )

  cpg_collections[["PANEL_APP"]] <- as.data.frame(
    gene_panels$records |>
      dplyr::left_join(
        dplyr::select(
          gene_info, entrezgene, gene_biotype
        ),
        by = "entrezgene"
      ) |>
      dplyr::select(entrezgene, gene_biotype, 
                    gepa_moi, gepa_phenotype) |>
      dplyr::distinct() |>
      dplyr::group_by(entrezgene, gene_biotype) |>
      dplyr::summarise(
        moi = paste(
          unique(sort(gepa_moi)), 
          collapse = "&"),
        phenotypes = paste(
          unique(sort(gepa_phenotype)),
          collapse = ";"
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(phenotypes = stringr::str_replace_all(
        phenotypes, "^(\\s{0,},\\s{1,})", ""
      )) |>
      dplyr::mutate(moi = stringr::str_replace(moi, "&NA$", "")) |>
      dplyr::mutate(source = "PANEL_APP") |>
      dplyr::mutate(moi = dplyr::if_else(
        stringr::str_detect(moi, "&|/"),
        "AD/AR",
        as.character(moi)
      )) |>
      dplyr::mutate(moi = dplyr::if_else(
        moi == "NA",
        as.character(NA),
        as.character(moi)
      ))
  )

  all_predisposition_incidental <- as.data.frame(
    #cpg_collections[["CGC"]] |>
    cpg_collections[['CANVAR_UK']] |>
      #dplyr::bind_rows(cpg_collections[['CANVAR_UK']]) |>
      dplyr::bind_rows(cpg_collections[["TCGA_PANCAN_2018"]]) |>
      dplyr::bind_rows(cpg_collections[["CURATED_OTHER"]]) |>
      dplyr::bind_rows(cpg_collections[["PANEL_APP"]]) |>
      dplyr::bind_rows(cpg_collections[["ACMG_SF"]]) |>
      dplyr::bind_rows(cpg_collections[["CPIC_PGX_ONCOLOGY"]]) |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::mutate(moi = dplyr::if_else(
        is.na(moi), "", as.character(moi)
      )) |>
      dplyr::mutate(predisp_syndrome_cui = dplyr::if_else(
        is.na(predisp_syndrome_cui), "",
        as.character(predisp_syndrome_cui)
      )) |>
      dplyr::mutate(susceptibility_cui = dplyr::if_else(
        is.na(susceptibility_cui), "",
        as.character(susceptibility_cui)
      )) |>
      dplyr::mutate(phenotypes = dplyr::if_else(
        is.na(phenotypes), "",
        as.character(phenotypes)
      )) |>
      dplyr::left_join(
        dplyr::select(
          gene_info,
          symbol, 
          entrezgene
        ),
        by = c("entrezgene"), multiple = "all"
      ) |>
      
      ## black list genes - no prominent evidence
      ## for cancer relevance
      ## 1) SERPINA1, GJB2, DHCR7, ASPM, UROD, SLC25A13 
      dplyr::filter(
        !(entrezgene %in% c(5265, 2706, 1717, 259266,
                            10165, 7389))
      ) |>
      
      dplyr::group_by(symbol, entrezgene, gene_biotype) |>
      dplyr::summarise(
        moi = paste(unique(sort(moi)), collapse = "&"),
        predisp_syndrome_cui = paste(
          unique(predisp_syndrome_cui),
          collapse = "&"
        ),
        susceptibility_cui = paste(
          unique(susceptibility_cui),
          collapse = "&"
        ),
        predisp_source = paste(unique(sort(source)), collapse = "&"),
        phenotypes = paste(unique(sort(phenotypes)), collapse = "; "),
        .groups = "drop"
      ) |>
      dplyr::mutate(phenotypes = stringr::str_replace_all(
        phenotypes, "(\\|NA)|(; NA)|(^;)", ""
      )) |>
      dplyr::mutate(phenotypes = stringr::str_replace_all(
        phenotypes, "^((( ){1,})?;?)", ""
      )) |>
      dplyr::mutate(moi = stringr::str_replace_all(moi, "^&|&$", "")) |>
      dplyr::mutate(predisp_syndrome_cui = stringr::str_replace_all(
        predisp_syndrome_cui, "^&|&$", ""
      )) |>
      dplyr::mutate(predisp_syndrome_cui = stringr::str_replace_all(
        predisp_syndrome_cui, ", ", "&"
      )) |>
      
      ## Make HFE ACMG_SF only
      dplyr::mutate(predisp_source = dplyr::if_else(
        entrezgene == 3077,
        "ACMG_SF",
        as.character(predisp_source)
      )) |>
      dplyr::mutate(susceptibility_cui = stringr::str_replace_all(
        susceptibility_cui, "^&|&$", ""
      )) |>
      dplyr::mutate(moi = dplyr::if_else(
        stringr::str_detect(moi, "AD/AR"),
        "AD/AR", as.character(moi)
      )) |>
      dplyr::mutate(moi = dplyr::if_else(
        !stringr::str_detect(moi, "AD/AR") &
          stringr::str_detect(moi, "AR") &
          stringr::str_detect(moi, "AD"),
        "AD/AR",
        as.character(moi)
      )) |>
      dplyr::left_join(
        mod_moi_predisposition, 
        by = "entrezgene", multiple = "all") |>
      dplyr::left_join(
        clingen_moi, 
        by = "entrezgene", multiple = "all") |>
      #dplyr::rename(mechanism_of_disease = mod_maxwell) |>
      dplyr::mutate(mechanism_of_disease = dplyr::case_when(
        !is.na(mod_maxwell) & nchar(mod_maxwell) > 0 ~ mod_maxwell,
        (is.na(mod_maxwell) | nchar(mod_maxwell) == 0) & !is.na(mod_clingen) ~ mod_clingen,
        TRUE ~ as.character(NA)
      )) |>
      dplyr::mutate(moi = dplyr::case_when(
        nchar(moi) == 0 & !is.na(moi_maxwell) ~ moi_maxwell,
        nchar(moi) == 0 & is.na(moi_clingen) & is.na(moi_maxwell) ~ as.character(NA),
        nchar(moi) == 0 & is.na(moi_maxwell) & !is.na(moi_clingen) ~ moi_clingen,
        TRUE ~ as.character(moi)
      )) |>
      dplyr::mutate(mechanism_of_disease = dplyr::case_when(
        ## manually curated - evidence for LoF as mechanism of disease
        .data$entrezgene == 8216 ~ "LoF", #LZTR1
          .data$entrezgene == 80169 ~ "LoF", #CTC1
          .data$entrezgene == 51750 ~ "LoF", #RTEL1
        TRUE ~ as.character(.data$mechanism_of_disease)
      )) |>
      dplyr::mutate(predisp_syndrome_cui = dplyr::if_else(
        nchar(predisp_syndrome_cui) == 0,
        as.character(NA),
        as.character(predisp_syndrome_cui)
      )) |>
      dplyr::mutate(susceptibility_cui = dplyr::if_else(
        nchar(susceptibility_cui) == 0,
        as.character(NA),
        as.character(susceptibility_cui)
      )) |>
      dplyr::mutate(moi = dplyr::if_else(
        nchar(moi) == 0, as.character(NA), as.character(moi)
      )) |>
      dplyr::rename(predisp_cancer_cui = susceptibility_cui) |>
      dplyr::select(-c("moi_maxwell","moi_clingen",
                       "mod_maxwell","mod_clingen",
                       "disease_label","disease_id_mondo",
                       "haploinsufficiency_description",
                       "haploinsufficiency_score")) |>
      dplyr::mutate(predisp_source = dplyr::if_else(
        stringr::str_detect(predisp_source, "CURATED_OTHER") &
          stringr::str_detect(predisp_source, "&"),
        stringr::str_replace(predisp_source, "^CURATED_OTHER&|&CURATED_OTHER", ""),
        as.character(predisp_source)
      )) |>
      dplyr::rename(
        cpg_source = predisp_source,
        cpg_syndrome_cui = predisp_syndrome_cui,
        cpg_cancer_cui = predisp_cancer_cui,
        cpg_mod = mechanism_of_disease,
        cpg_moi = moi,
        cpg_phenotypes = phenotypes
      ) |>
      dplyr::mutate(
        cpg_phenotypes = stringr::str_replace_all(
          cpg_phenotypes, 
          "[\r\n\t]", 
          "")
      )
  )
  
  return(all_predisposition_incidental)
}
