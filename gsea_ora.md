------------------------------------------------------------------------

author: “Aaron Mohammed” output: github_document —

``` r
library(readr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(simplifyEnrichment)
library(openxlsx)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(devEMF)

getOption("clusterProfiler.download.method")
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

project_dir <- getwd()

diff_path <- file.path(project_dir,
                       "DGE_output",
                       "tables",
                       "differential")

abundance_path <- file.path(project_dir,
                            "DGE_output",
                            "tables",
                            "processed_abundance")

salmon_counts_df <- read_tsv(file.path(project_dir,
                             "rnaseq_output",
                             "star_salmon",
                             "salmon.merged.gene_counts.tsv"))

annots_df <- salmon_counts_df[ ,c("gene_id", "gene_name")]
colnames(annots_df) <- c("ENSEMBL_ID.version", "GENE_NAME")
annots_df <- annots_df[,c("GENE_NAME", "ENSEMBL_ID.version")]

annots_df$ENSEMBL_ID <- sub("(.*?)(\\..*)?$", "\\1", annots_df$ENSEMBL_ID.version)
  
annots_df$ENTREZ_ID <- mapIds(x = org.Mm.eg.db,
                keys = annots_df$GENE_NAME,
                column = "ENTREZID",
                keytype = "SYMBOL",
                multiVals = "first")

temp = list.files(path= file.path(diff_path), 
                  pattern="\\.results.tsv$", full.names = TRUE)
deseq2_files = lapply(temp, read_tsv)

contrast_names <- sub("\\.deseq2.results.tsv.*", "", 
                      list.files(path= file.path(diff_path), pattern="\\.results.tsv$"))

names(deseq2_files) <- contrast_names
```

``` r
for (i in 1:length(deseq2_files)) {
  gene_ids <- deseq2_files[[i]]$gene_id
  gene_names <- annots_df[which(annots_df$ENSEMBL_ID.version %in% gene_ids), "GENE_NAME"]
  cnames <- colnames(deseq2_files[[i]])
  deseq2_files[[i]][,"gene_name"] <- gene_names
  deseq2_files[[i]] <- deseq2_files[[i]][,c("gene_name",cnames)]
}
```

``` r
expression <- read_tsv(file.path(abundance_path, "all.vst.tsv"))
cnames <- colnames(expression)

expression[,"gene_name"] <- annots_df[which(annots_df$ENSEMBL_ID.version %in% expression$gene_id),
                                      "GENE_NAME"]

expression <- expression[,c("gene_name", cnames)]

# write.xlsx(expression,
#            file.path(project_dir,"expression_vst.xlsx"),
#            rowNames = FALSE)
```

``` r
ensembl.to.symbol <- function (enrich_list, annots_df) {
  for (i in 1:length(enrich_list)) {
    for (j in 1:length(enrich_list[[i]][["core_enrichment"]])) {
      eids <- enrich_list[[i]][["core_enrichment"]][j]
      eids2 <- enrich_list[[i]][["Significant_DEGs"]][j]
      
      # Split the vector into individual IDs
      ids <- strsplit(eids, "/")[[1]]
      ids2 <- strsplit(eids2, "/")[[1]]
      
      # Match and replace IDs with gene names
      gene_names <- sapply(ids, function(id) {
        gene_name <- annots_df$GENE_NAME[match(id, annots_df$ENSEMBL_ID)]
        if (is.na(gene_name)) id else gene_name
      })
      
      gene_names2 <- sapply(ids2, function(id) {
        gene_name <- annots_df$GENE_NAME[match(id, annots_df$ENSEMBL_ID)]
        if (is.na(gene_name)) id else gene_name
      })
      
      # Recombine gene names back into a single string
      gene_names <- paste(gene_names, collapse = "/")
      gene_names -> enrich_list[[i]][["core_enrichment"]][j]
      
      gene_names2 <- paste(gene_names2, collapse = "/")
      gene_names2 -> enrich_list[[i]][["Significant_DEGs"]][j]
      }
    }
  return(enrich_list)
}
```

``` r
inter.func <- function (enri_obj, sig_genes) {
  splited <- strsplit(enri_obj@result$core_enrichment, split = "/")
  names(splited) <- enri_obj@result$ID
  
  intersected <- {}
  for (i in 1:length(splited)) {
    intersected[[i]] <- splited[[i]][which(splited[[i]] %in% sig_genes)]
  }
  
  names(intersected) <- enri_obj@result$ID
  
  concatenated <- vector(mode = "list", length = length(intersected))
  for (i in 1:length(intersected)) {
    concatenated[i] <- paste(intersected[[i]], collapse = "/")
  }
  
  names(concatenated) <- enri_obj@result$ID
  
  enri_obj@result$Significant_DEGs <- unlist(concatenated)
  
  return(enri_obj)
}
```

``` r
go.gsea <- function(ranked_gene_list, up_dir, down_dir, groupName, db, 
                    padj_cutoff, color_by, minGSS, annots_df, sig_genes) {
  
  set.seed(123)
  
  up_gsea_go <- file.path(up_dir, "GO")
  dir.create(up_gsea_go)
  
  down_gsea_go <- file.path(down_dir, "GO")
  dir.create(down_gsea_go)
  
  onts <- c("BP", "CC", "MF")
  
  exc_BP <- FALSE
  exc_CC <- FALSE
  exc_MF <- FALSE
  
  # GSEA GO
  gsea_GO_BP <- gseGO(geneList=ranked_gene_list,
                      ont ="BP", 
                      keyType = "ENSEMBL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = padj_cutoff, # cutoff for adjusted pvalue
                      minGSSize = minGSS)
  
  if (nrow(gsea_GO_BP@result) == 0) {
    exc_BP <- TRUE
  }

  gsea_GO_CC <- gseGO(geneList=ranked_gene_list,
                      ont ="CC", 
                      keyType = "ENSEMBL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = padj_cutoff, # cutoff for adjusted pvalue
                      minGSSize = minGSS)
  
  if (nrow(gsea_GO_CC@result) == 0) {
    exc_CC <- TRUE
  }
  
  gsea_GO_MF <- gseGO(geneList=ranked_gene_list,
                      ont ="MF", 
                      keyType = "ENSEMBL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = padj_cutoff, # cutoff for adjusted pvalue
                      minGSSize = minGSS)
  
  if (nrow(gsea_GO_MF@result) == 0) {
    exc_MF <- TRUE
  }
  
  exc <- c(exc_BP, exc_CC, exc_MF)
  gse_results <- list(gsea_GO_BP, gsea_GO_CC, gsea_GO_MF)
  names(gse_results) <- onts
  
  gse_results <- gse_results[which(!exc)]
  onts <- onts[which(!exc)]
  
    if (length(gse_results) > 1) {
      
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                  # UP
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
      for (x in 1:length(gse_results)) {
          gse_results[[x]] <- inter.func(gse_results[[x]], sig_genes)
      }
      
      up_gsea_GO_list <- vector(mode = "list", length = length(gse_results))
      names(up_gsea_GO_list) <- onts
      
      for (x in onts) {
        up_gsea_GO_list[[x]] <- subset_enrichResult(gse_results[[x]],
                                                    which(gse_results[[x]]@result$NES > 0))
      }
      
      for (x in names(up_gsea_GO_list)) {
        if (nrow(up_gsea_GO_list[[x]]@result) == 0) {
          up_gsea_GO_list[[x]] <- NULL
        }
      }
      
        if (length(up_gsea_GO_list) > 1) {
          UP_list <- vector(mode = "list", length = length(up_gsea_GO_list))
          names(UP_list) <- names(up_gsea_GO_list)
          
          for (x in names(up_gsea_GO_list)) {
            UP_list[[x]] <- up_gsea_GO_list[[x]]@result
          }
          write.xlsx(UP_list, file.path(up_gsea_go, 
                                        paste0("UP_GSEA_GO", sep = "_",
                                               groupName, sep="_", 
                                               "ENSEMBL_IDs.xlsx")))
          
          UP_list <- ensembl.to.symbol(UP_list, annots_df)
          
          write.xlsx(UP_list, file.path(up_gsea_go, paste0("UP_GSEA_GO", sep = "_",
                                                           groupName, sep="_", 
                                                           "Gene_Names.xlsx")))
          
          for (o in names(up_gsea_GO_list)) {
            if (o == "BP") {
              title <- "Biological Processes"
              } else if (o == "MF") {
              title <- "Molecular Functions"
              } else if (o == "CC") {
              title <- "Cellular Components"
              }
            if (nrow(up_gsea_GO_list[[o]]@result) >= 20) {
              h <- 9
              if (max(nchar(up_gsea_GO_list[[o]]@result$Description[1:20])) >= 120) {
                h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description[1:20]))/45)}
              } else if (nrow(up_gsea_GO_list[[o]]@result) >= 10 & 
                         nrow(up_gsea_GO_list[[o]]@result) < 20) {
                h <- 7
              if (max(nchar(up_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description))/45)}} 
            else if (nrow(up_gsea_GO_list[[o]]@result) >= 5 & 
                     nrow(up_gsea_GO_list[[o]]@result) < 10) {
              h <- 5
              if (max(nchar(up_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description))/45)}} 
            else {
              h <- 4
              if (max(nchar(up_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description))/45)}}
            
            
            dp <- dotplot(up_gsea_GO_list[[o]], 
                          showCategory=20, 
                          label_format = 60,
                          color = color_by) + 
              ggtitle(title) +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            
            
            pdf(height = h,
                width = 12,
                file.path(up_gsea_go, 
                          paste0("UP", sep="_",
                                 o, sep="_",
                                 groupName, sep="_", 
                                 "GSEA_GO_dotplot.pdf")))
            print(dp)
            dev.off()

            set.seed(123)
            
            pt <- pairwise_termsim(up_gsea_GO_list[[o]])
            
            emap <- emapplot(pt, 
                             color = color_by) + 
              ggtitle(title) +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
            pdf(height = 15,
                width=15,
                file.path(up_gsea_go, 
                          paste0("UP", sep="_",
                                 o, sep="_",
                                 groupName, sep="_", 
                                 "GSEA_GO_emap.pdf")))
            print(emap)
            dev.off()

          }
        }
      
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                  # DOWN
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
        
      down_gsea_GO_list <- vector(mode = "list", length = length(gse_results))
      names(down_gsea_GO_list) <- onts
      
      for (x in onts) {
        down_gsea_GO_list[[x]] <- subset_enrichResult(gse_results[[x]],
                                                    which(gse_results[[x]]@result$NES < 0))
      }
      
      for (x in names(down_gsea_GO_list)){
        if (nrow(down_gsea_GO_list[[x]]@result) == 0) {
          down_gsea_GO_list[[x]] <- NULL
        }
      }
      
        if (length(down_gsea_GO_list) > 1) {
          DOWN_list <- vector(mode = "list", length = length(down_gsea_GO_list))
          names(DOWN_list) <- names(down_gsea_GO_list)
          
          for (x in names(down_gsea_GO_list)) {
            DOWN_list[[x]] <- down_gsea_GO_list[[x]]@result
          }
          write.xlsx(DOWN_list, file.path(down_gsea_go, 
                                        paste0("DOWN_GSEA_GO", sep = "_",
                                               groupName, sep="_", 
                                               "ENSEMBL_IDs.xlsx")))
          
          DOWN_list <- ensembl.to.symbol(DOWN_list, annots_df)
          
          write.xlsx(DOWN_list, file.path(down_gsea_go, paste0("DOWN_GSEA_GO", sep = "_",
                                                           groupName, sep="_", 
                                                           "Gene_Names.xlsx")))
          
          for (o in names(down_gsea_GO_list)) {
            if (o == "BP") {
              title <- "Biological Processes"
            } else if (o == "MF") {
              title <- "Molecular Functions"
            } else if (o == "CC") {
              title <- "Cellular Components"
            }
            
            if (nrow(down_gsea_GO_list[[o]]@result) >= 20) {
              h <- 9
              if (max(nchar(down_gsea_GO_list[[o]]@result$Description[1:20])) >= 120) {
                h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description[1:20]))/45)
              }
            } else if (nrow(down_gsea_GO_list[[o]]@result) >= 10 & 
                       nrow(down_gsea_GO_list[[o]]@result) < 20) {
              h <- 7
              if (max(nchar(down_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description))/45)
              }
            } else if (nrow(down_gsea_GO_list[[o]]@result) >= 5 & 
                       nrow(down_gsea_GO_list[[o]]@result) < 10) {
              h <- 5
              if (max(nchar(down_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description))/45)
              }
            } else {
              h <- 4
              if (max(nchar(down_gsea_GO_list[[o]]@result$Description)) >= 120) {
                h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description))/45)
              }
            }
            
            
            dp <- dotplot(down_gsea_GO_list[[o]], 
                          showCategory=20, 
                          label_format = 60,
                          color = color_by) + 
              ggtitle(title) +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            
            
            pdf(height = h,
                width = 12,
                file.path(down_gsea_go, 
                          paste0("DOWN", sep="_",
                                 o, sep="_",
                                 groupName, sep="_", 
                                 "GSEA_GO_dotplot.pdf")))
            print(dp)
            dev.off()

            set.seed(123)
            
            pt <- pairwise_termsim(down_gsea_GO_list[[o]])
            
            emap <- emapplot(pt, 
                             color = color_by) + 
              ggtitle(title) +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
            pdf(height = 15,
                width=15,
                file.path(down_gsea_go, 
                          paste0("DOWN", sep="_",
                                 o, sep="_",
                                 groupName, sep="_", 
                                 "GSEA_GO_emap.pdf")))
            print(emap)
            dev.off()

          }
          
        }
    }
  
  length(list.files(up_gsea_go)) -> nfiles_up
  length(list.files(down_gsea_go)) -> nfiles_down
  
  if (nfiles_up == 0) {
    new_up_gsea_go <- paste0(up_gsea_go, "-none_enriched")
    file.rename(up_gsea_go, new_up_gsea_go)
  } else if (nfiles_up == 2) {
    new_up_gsea_go <- paste0(up_gsea_go, "-padj_thresh_not_met")
    file.rename(up_gsea_go, new_up_gsea_go)
  }
  
  if (nfiles_down == 0) {
    new_down_gsea_go <- paste0(down_gsea_go, "-none_enriched")
    file.rename(down_gsea_go, new_down_gsea_go)
  } else if (nfiles_down == 2) {
    new_down_gsea_go <- paste0(down_gsea_go, "-padj_thresh_not_met")
    file.rename(down_gsea_go, new_down_gsea_go)
  }
  
}
```

``` r
entrez.to.symbol <- function (KEGG_result, annots_df) {
  for (i in 1:length(KEGG_result$core_enrichment)) {
      eids <- KEGG_result$core_enrichment[i]
      eids2 <- KEGG_result$Significant_DEGs
      
      # Split the vector into individual IDs
      ids <- strsplit(eids, "/")[[1]]
      ids2 <- strsplit(eids2, "/")[[1]]
      
      # Match and replace IDs with gene names
      gene_names <- sapply(ids, function(id) {
        gene_name <- annots_df$GENE_NAME[match(id, annots_df$ENTREZ_ID)]
        if (is.na(gene_name)) id else gene_name
      })
      
      gene_names2 <- sapply(ids2, function(id) {
        gene_name <- annots_df$GENE_NAME[match(id, annots_df$ENTREZ_ID)]
        if (is.na(gene_name)) id else gene_name
      })
      
      # Recombine gene names back into a single string
      gene_names <- paste(gene_names, collapse = "/")
      gene_names -> KEGG_result$core_enrichment[i]
      
      gene_names2 <- paste(gene_names2, collapse = "/")
      gene_names2 -> KEGG_result$Significant_DEGs
      
    }
  return(KEGG_result)
}
```

``` r
kegg.gsea <- function(ranked_gene_list, up_dir, down_dir, groupName, db, organism, 
                      padj_cutoff, color_by, minGSS, annots_df, sig_genes) {
  
  set.seed(123)
  
  up_gsea_kegg <- file.path(up_dir, "KEGG")
  dir.create(up_gsea_kegg)
  
  down_gsea_kegg <- file.path(down_dir, "KEGG")
  dir.create(down_gsea_kegg)
  
  gsea_KEGG <- gseKEGG(geneList=ranked_gene_list,
                       organism = organism,
                       use_internal_data = F,
                       pAdjustMethod = "BH",
                       pvalueCutoff = padj_cutoff, # cutoff for adjusted pvalue
                       minGSSize = minGSS)
  
  if (nrow(gsea_KEGG@result) > 0) {
    
    gsea_KEGG <- inter.func(gsea_KEGG, sig_genes)
    
    up_gsea_KEGG <- subset_enrichResult(gsea_KEGG, 
                                        which(gsea_KEGG@result$NES > 0))
    
    if (nrow(up_gsea_KEGG@result) > 0) {
      
      write.xlsx(up_gsea_KEGG@result, 
                 file.path(up_gsea_kegg, 
                           paste0("UP_GSEA_KEGG", sep = "_",
                                  groupName, sep="_", 
                                  "ENTREZ_IDs.xlsx")))
      
      up_gsea_KEGG@result <- entrez.to.symbol(up_gsea_KEGG@result, annots_df)
      
      write.xlsx(up_gsea_KEGG@result, 
                 file.path(up_gsea_kegg, 
                           paste0("UP_GSEA_KEGG", sep = "_",
                                  groupName, sep="_", 
                                  "Gene_Names.xlsx")))
      
    if (nrow(up_gsea_KEGG@result) >= 20) {
      h <- 9
      if (max(nchar(up_gsea_KEGG@result$Description[1:20])) >= 120) {
        h <- h + floor(max(nchar(up_gsea_KEGG@result$Description[1:10]))/45)
      }
    } else if (nrow(up_gsea_KEGG@result) >= 10 & nrow(up_gsea_KEGG@result) < 20) {
      h <- 7
      if (max(nchar(up_gsea_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(up_gsea_KEGG@result$Description))/45)
      }
    } else if (nrow(up_gsea_KEGG@result) >= 5 & nrow(up_gsea_KEGG@result) < 10) {
      h <- 5
      if (max(nchar(up_gsea_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(up_gsea_KEGG@result$Description))/45)
      }
    } else {
      h <- 4
      if (max(nchar(up_gsea_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(up_gsea_KEGG@result$Description))/45)
      }
    }
      
      dp <- dotplot(up_gsea_KEGG, 
                    showCategory=20, 
                    label_format = 60,
                    color = color_by)
      
      pdf(height = h, 
          width = 12,
          file.path(up_gsea_kegg, 
                    paste0("UP", sep="_",
                           groupName, sep="_", 
                           "GSEA_KEGG_dotplot.pdf")))
      print(dp)
      dev.off()

      set.seed(123)
      
      pt <- pairwise_termsim(up_gsea_KEGG)
      emap <- emapplot(pt, 
                       color = color_by)
      
      pdf(height = 15,
          width=15,
          file.path(up_gsea_kegg, 
                    paste0("UP", sep="_",
                           groupName, sep="_", 
                           "GSEA_KEGG_emap.pdf")))
      print(emap)
      dev.off()

    }
    
    down_gsea_KEGG <- subset_enrichResult(gsea_KEGG, 
                                          which(gsea_KEGG@result$NES < 0))
    
    if (nrow(down_gsea_KEGG@result) > 0) {
      
      write.xlsx(down_gsea_KEGG@result, 
                 file.path(down_gsea_kegg, 
                           paste0("DOWN_GSEA_KEGG", sep = "_",
                                  groupName, sep="_", 
                                  "ENTREZ_IDs.xlsx")))
      
      down_gsea_KEGG@result <- entrez.to.symbol(down_gsea_KEGG@result, annots_df)
      
      write.xlsx(down_gsea_KEGG@result, 
                 file.path(down_gsea_kegg, 
                           paste0("DOWN_GSEA_KEGG", sep = "_",
                                  groupName, sep="_", 
                                  "Gene_Names.xlsx")))
      
      if (nrow(down_gsea_KEGG@result) >= 20) {
        h <- 9
        if (max(nchar(down_gsea_KEGG@result$Description[1:20])) >= 120) {
          h <- h + floor(max(nchar(down_gsea_KEGG@result$Description[1:10]))/45)
        }
      } else if (nrow(down_gsea_KEGG@result) >= 10 & nrow(down_gsea_KEGG@result) < 20) {
        h <- 7
        if (max(nchar(down_gsea_KEGG@result$Description)) >= 120) {
          h <- h + floor(max(nchar(down_gsea_KEGG@result$Description))/45)
        }
      } else if (nrow(down_gsea_KEGG@result) >= 5 & nrow(down_gsea_KEGG@result) < 10) {
        h <- 5
        if (max(nchar(down_gsea_KEGG@result$Description)) >= 120) {
          h <- h + floor(max(nchar(down_gsea_KEGG@result$Description))/45)
        }
      } else {
        h <- 4
        if (max(nchar(down_gsea_KEGG@result$Description)) >= 120) {
          h <- h + floor(max(nchar(down_gsea_KEGG@result$Description))/45)
        }
      }
      
      dp <- dotplot(down_gsea_KEGG, 
                    showCategory=20, 
                    label_format = 60,
                    color = color_by)
      
      pdf(height = h, 
          width=12,
          file.path(down_gsea_kegg, 
                    paste0("DOWN", sep="_",
                           groupName, sep="_", 
                           "GSEA_KEGG_dotplot.pdf")))
      print(dp)
      dev.off()

      set.seed(123)
      
      pt <- pairwise_termsim(down_gsea_KEGG)
      
      emap <- emapplot(pt, 
                       color = color_by)
      
      pdf(height = 15, 
          width=15,
          file.path(down_gsea_kegg, 
                    paste0("DOWN", sep="_",
                           groupName, sep="_", 
                           "GSEA_KEGG_emap.pdf")))
      print(emap)
      dev.off()

    }
  }
  
  length(list.files(up_gsea_kegg)) -> nfiles_up
  length(list.files(down_gsea_kegg)) -> nfiles_down
  
  if (nfiles_up == 0) {
    new_up_gsea_kegg <- paste0(up_gsea_kegg, "-none_enriched")
    file.rename(up_gsea_kegg, new_up_gsea_kegg)
  } else if (nfiles_up == 2) {
    new_up_gsea_kegg <- paste0(up_gsea_kegg, "-padj_thresh_not_met")
    file.rename(up_gsea_kegg, new_up_gsea_kegg)
  }
  
  if (nfiles_down == 0) {
    new_down_gsea_kegg <- paste0(down_gsea_kegg, "-none_enriched")
    file.rename(down_gsea_kegg, new_down_gsea_kegg)
  } else if (nfiles_down == 2) {
    new_down_gsea_kegg <- paste0(down_gsea_kegg, "-padj_thresh_not_met")
    file.rename(down_gsea_kegg, new_down_gsea_kegg)
  }
  
}
```

``` r
go.ora.plots <- function (enrichResult, ont, path, up_or_down, groupName, 
                          color_by, padj_cutoff) {
  
  ont_sub <- subset_enrichResult(enrichResult, 
                                 which(enrichResult@result$ONTOLOGY==ont & 
                                         enrichResult@result$p.adjust <= padj_cutoff))
  
  if (nrow(ont_sub@result) > 0) {
  
  if (ont == "BP") {
    title <- "Biological Processes"
    } else if (ont == "MF") {
      title <- "Molecular Functions"
      } else if (ont == "CC") {
        title <- "Cellular Components"
        }
  
   if (nrow(enrichResult@result) >= 20) {
     h <- 8
     if (max(nchar(enrichResult@result$Description[1:20])) >= 120) {
       h <- h + floor(max(nchar(enrichResult@result$Description[1:20]))/45)
       }
     } else if (nrow(enrichResult@result) >= 10 & 
                nrow(enrichResult@result) < 20) {
       h <- 6
       if (max(nchar(enrichResult@result$Description)) >= 120) {
        h <- h + floor(max(nchar(enrichResult@result$Description))/45)
      }
    } else if (nrow(enrichResult@result) >= 5 & nrow(enrichResult@result) < 10) {
      h <- 5
      if (max(nchar(enrichResult@result$Description)) >= 120) {
        h <- h + floor(max(nchar(enrichResult@result$Description))/45) 
      }
    } else {
      h <- 4
      if (max(nchar(enrichResult@result$Description)) >= 120) {
        h <- h + floor(max(nchar(enrichResult@result$Description))/45)
      }
    }
  
  bp <- barplot(ont_sub,
        showCategory = 20,
        color = color_by,
        label_format = 60) +
  ggtitle(title) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  pdf(height = h, 
      width=12,
      file.path(path, paste0(up_or_down, sep="_", 
                            groupName, sep="_",
                            "ORA", sep= "_",
                            "GO", sep= "_",
                            ont, sep= "_",
                            "barplot.pdf")))
  
  print(bp)
  dev.off()

  set.seed(123)
  
  pt <- pairwise_termsim(ont_sub)
  
  emap <- emapplot(pt, color = color_by) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = "bold"))
  
  pdf(height = 15, 
      width=15,
      file.path(path, paste0(up_or_down, sep="_", 
                            groupName, sep="_",
                            "ORA", sep= "_",
                            "GO", sep= "_",
                            ont, sep= "_",
                            "emap.pdf")))
  print(emap)
  dev.off()

  }
  
}
```

``` r
go.ora <- function(df, path, up_or_down, groupName, db, color_by, padj_cutoff, minGSS) {
  
  set.seed(123)
  
  ora_go_path <- file.path(path, "GO")
  dir.create(ora_go_path)


  # Get names of the genes
  gene_list <- df$gene_name
  
  # ORA GO
  ora_GO <- enrichGO(gene = gene_list,
                     OrgDb= db,
                     keyType = "SYMBOL",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = padj_cutoff, # cutoff for adjusted pvalue
                     minGSSize = minGSS) 
  
    if (is.null(ora_GO)) {
      skip <- TRUE
      } else if (nrow(ora_GO@result) == 0) {
        skip <- TRUE
      } else if (min(ora_GO@result$p.adjust) > padj_cutoff) {
        skip <- TRUE
      } else {skip <- FALSE}
  
  if (!skip) {
    onts <- unique(ora_GO@result$ONTOLOGY)
    
      ora_GO_list <- vector(mode = "list", length = length(onts))
      names(ora_GO_list) <- onts
      
  for (o in onts) {
    ora_GO_list[[o]] <- ora_GO@result[which(ora_GO@result$ONTOLOGY == o), ]
  }
  
  write.xlsx(ora_GO_list, 
             file.path(ora_go_path, paste0(up_or_down, sep="_", 
                                   groupName, sep="_", 
                                   "ORA_GO.xlsx")))
   for (o in onts) {
    
    go.ora.plots(enrichResult = ora_GO, 
              ont = o, 
              path = ora_go_path, 
              up_or_down = up_or_down, 
              groupName = groupName,
              color_by = color_by,
              padj_cutoff = padj_cutoff)
  }
  
  }
  
    length(list.files(ora_go_path)) -> nfiles
  
    if (nfiles == 0) {
    path_new <- paste0(ora_go_path, "-none_enriched")
    file.rename(ora_go_path, path_new)
  } else if (nfiles == 1) {
    path_new <- paste0(ora_go_path, "-padj_thresh_not_met")
    file.rename(ora_go_path, path_new)
  }
  
}
```

``` r
kegg.ora <- function(df, up_or_down, path, groupName, db, organism, color_by, 
                     padj_cutoff, minGSS) {
  
  set.seed(123)
  
  ora_kegg_path <- file.path(path, "KEGG")
  dir.create(ora_kegg_path, 
             showWarnings = FALSE)

  # Get names of the genes
  gene_list <- df$gene_name
  
  # Convert the names of the genes to ENTREZID for KEGG
  entrz_id <- mapIds(x = db,
                     keys = df$gene_id,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
  
  names(entrz_id) <- NULL
  entrz_id <- entrz_id[-which(is.na(entrz_id))]
  
  # ORA KEGG
  ora_KEGG <- enrichKEGG(gene = entrz_id,
                 organism = organism,
                 use_internal_data = F,
                 pAdjustMethod = "BH",
                 pvalueCutoff = padj_cutoff,
                 minGSSize = minGSS)
  
    if (is.null(ora_KEGG)) {
      skip <- TRUE
      } else if (nrow(ora_KEGG@result) == 0) {
        skip <- TRUE
      } else if (min(ora_KEGG@result$p.adjust) > padj_cutoff) {
        skip <- TRUE
      } else {skip <- FALSE}
  
  if (!skip) {
  
  write.xlsx(ora_KEGG@result, file.path(ora_kegg_path, paste0(up_or_down, sep="_", 
                                                   groupName, sep="_", 
                                                   "ORA_KEGG.xlsx")))
    
   if (nrow(ora_KEGG@result) >= 20) {
     h <- 8
     if (max(nchar(ora_KEGG@result$Description[1:20])) >= 120) {
       h <- h + floor(max(nchar(ora_KEGG@result$Description[1:20]))/45)}
     }
    
    if (nrow(ora_KEGG@result) >= 10 & nrow(ora_KEGG@result) < 20) {
      h <- 6
      if (max(nchar(ora_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(ora_KEGG@result$Description))/45)
      }
    }
    
    if (nrow(ora_KEGG@result) >= 5 & nrow(ora_KEGG@result) < 10) {
      h <- 5
      if (max(nchar(ora_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(ora_KEGG@result$Description))/45)
        } else {
          h <- 4
          if (max(nchar(ora_KEGG@result$Description)) >= 120) {
            h <- h + floor(max(nchar(ora_KEGG@result$Description))/45)}
        }
    }
  
    bp <- barplot(ora_KEGG, 
                  showCategory = 20,
                  color = color_by,
                  label_format = 60)
    
    pdf(file.path(ora_kegg_path, 
                  paste0(up_or_down, sep="_", groupName, sep="_", "ORA_KEGG_barplot.pdf")), 
        height = h, width=10)
    print(bp)
    dev.off()

    pt <- pairwise_termsim(ora_KEGG)
    
    emap <- emapplot(pt, color = color_by)
    
    pdf(file.path(ora_kegg_path, 
                  paste0(up_or_down, sep="_", groupName, sep="_", "ORA_KEGG_emap.pdf")), 
        height = 15, 
        width=15)
    print(emap)
    dev.off()
    
  }
  
  
    length(list.files(ora_kegg_path)) -> nfiles
  
    if (nfiles == 0) {
    new_ora_kegg_path <- paste0(ora_kegg_path, "-none_enriched")
    file.rename(ora_kegg_path, new_ora_kegg_path)
  } else if (nfiles == 1) {
    new_ora_kegg_path <- paste0(ora_kegg_path, "-padj_thresh_not_met")
    file.rename(ora_kegg_path, new_ora_kegg_path)
  }

}
```

``` r
gsea.ora <- function (deseq2_files, proj_path, dir_name = "GSEA_ORA", log2FC_cutoff, 
                      padj_cutoff, color_by, db, organism, minGSS, 
                      annots_df, expression, sample_groups, top_DGEs = 50) {
  
  path <- file.path(proj_path, "GSEA_ORA")
  dir.create(path, showWarnings = FALSE)

  for (comp in names(deseq2_files)) {
    
    result_dir <- file.path(path, comp)
    dir.create(result_dir)
    
    gsea_dir <- file.path(result_dir, "GSEA")
    dir.create(gsea_dir)

    ora_dir <- file.path(result_dir, "ORA")
    dir.create(ora_dir)

    df <- deseq2_files[[comp]]
    df$gene_id <- sub("(.*?)(\\..*)?$", "\\1", df$gene_id)
    
    groupName <- comp
    
    write.xlsx(df,
           file.path(file.path(result_dir),
                     paste0(groupName, sep="_",
                            "deseq2_result", sep="",
                            ".xlsx")),
           rowNames = FALSE)

    df <- df[order(-df$log2FoldChange),]

    write.xlsx(df,
           file.path(result_dir,
                     paste0(groupName, sep="_",
                            "deseq2_result", sep="_",
                            "sorted_by_log2FoldChange", sep="",
                            ".xlsx")),
           rowNames = FALSE)

    #######################################################################
                              # Volcano Plots
    #######################################################################

    upper <- round(max(df$log2FoldChange)) + 1
    lower <- round(min(df$log2FoldChange)) + 1
    lim <- max(upper, abs(lower))    
      
    sub_title_padj <- paste0("log2FC cutoff =", sep = " " ,
                        log2FC_cutoff, sep= "     ",
                        "padj cutoff =", sep=" ",
                        padj_cutoff)
    
    plot <- EnhancedVolcano(df,
                            x="log2FoldChange", y="padj",
                            lab= df$gene_name,
                            pCutoff = padj_cutoff,
                            FCcutoff = log2FC_cutoff,
                            xlim = c(-lim,lim),
                            title = groupName,
                            subtitle = sub_title_padj,
                            legendLabels  = c("NS", 
                                              expression(Log[2] ~ FC), 
                                              "padj", 
                                              expression(padj ~ and ~ log[2] ~ FC)))
    
    pdf(file.path(result_dir, 
                  paste0(groupName, sep="_",
                         "vol_plot_padj", sep = "",
                         ".pdf")),
        height = 10,
        width = 12)
    print(plot)
    dev.off()

    #######################################################################
                         # Gene Set Enrichment Analysis
    #######################################################################
    
    up_gsea_dir <- file.path(gsea_dir, "UP")
    dir.create(up_gsea_dir)
    
    down_gsea_dir <- file.path(gsea_dir, "DOWN")
    dir.create(down_gsea_dir)

    gene_list <- df$log2FoldChange
    names(gene_list) <- sub("(.*?)(\\..*)?$", "\\1", df$gene_id)
  
    sig_genes <- df[which(df$log2FoldChange >= log2FC_cutoff &
                              df$padj <= padj_cutoff),]
  
    sig_genes <- rbind(sig_genes,
                         df[which(df$log2FoldChange <= log2FC_cutoff &
                              df$padj <= padj_cutoff),])
  
    sig_genes <- sig_genes$gene_id
  
    
    go.gsea(ranked_gene_list = gene_list,
            up_dir = up_gsea_dir,
            down_dir = down_gsea_dir,
            groupName = groupName,
            db = db,
            padj_cutoff = padj_cutoff,
            color_by = color_by,
            minGSS = minGSS,
            annots_df = annots_df,
            sig_genes = sig_genes)
    
    names(gene_list) <- mapIds(x = db,
                               keys = df$gene_id,
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")
  
      gene_list <- gene_list[-which(is.na(names(gene_list)))]
      gene_list <- sort(gene_list, decreasing = TRUE)
  
      sig_genes <- annots_df[which(annots_df$ENSEMBL_ID %in% sig_genes), "ENTREZ_ID"]
      sig_genes <- sig_genes$ENTREZ_ID[-which(is.na(sig_genes$ENTREZ_ID))]
      names(sig_genes) <- NULL
  
      kegg.gsea(ranked_gene_list = gene_list,
                up_dir = up_gsea_dir,
                down_dir = down_gsea_dir,
                groupName = groupName,
                db = db,
                padj_cutoff = padj_cutoff,
                color_by = "p.adjust",
                minGSS = minGSS,
                organism = organism,
                annots_df = annots_df,
                sig_genes = sig_genes)

      ###################################################################
                        # Over-Representation Analysis
      ###################################################################
      
      up_ora_dir <- file.path(ora_dir, "UP")
      dir.create(up_ora_dir)
      
      down_ora_dir <- file.path(ora_dir, "DOWN")
      dir.create(down_ora_dir)
      
      df_up <- df[which(df$log2FoldChange >= log2FC_cutoff & df$padj <= padj_cutoff),]
      df_up <- df_up[order(-df_up$log2FoldChange), ]
      
      write.xlsx(df_up,
           file.path(result_dir,
                     paste0("UP", sep = "_",
                            groupName, sep="",
                            ".xlsx")),
           rowNames = FALSE)
      

      
      df_down <- df[which(df$log2FoldChange <= -log2FC_cutoff & df$padj <= padj_cutoff),]
      df_down <- df_down[order(df_down$log2FoldChange), ]
      
      write.xlsx(df_down,
                 file.path(result_dir,
                           paste0("DOWN", sep = "_",
                                  groupName, sep="",
                                  ".xlsx")),
                 rowNames = FALSE)
      
      go.ora(df = df_up,
             up_or_down = "UP",
             path = up_ora_dir,
             groupName = groupName,
             db = db,
             padj_cutoff = padj_cutoff,
             color_by = color_by,
             minGSS = minGSS)

      kegg.ora(df = df_up,
             up_or_down = "UP",
             path = up_ora_dir,
             groupName = groupName,
             db = db,
             organism = organism,
             padj_cutoff = padj_cutoff,
             color_by = color_by,
             minGSS = minGSS)
      
    go.ora(df = df_down,
           up_or_down = "DOWN",
           path = down_ora_dir,
           groupName = groupName,
           db = db,
           padj_cutoff = padj_cutoff,
           color_by = color_by,
           minGSS = minGSS)

    kegg.ora(df = df_down,
             up_or_down = "DOWN",
             path = down_ora_dir,
             groupName = groupName,
             db = db,
             organism = organism,
             padj_cutoff = padj_cutoff,
             color_by = color_by,
             minGSS = minGSS)

  ##############################################################################

  hm_colors <- colorRampPalette(c("#4ba9e7", "#0F5B8D", "#043250", "black",
                                  "#741C03", "#A92F08", "#FC4911"))

  # pat <- gsub("_vs_control", "", groupName)
  # pat <- gsub("_", ".", pat)

  up_genes <- df_up$gene_name
  up_genes <- up_genes[1:top_DGEs]
  up_exp <- expression[match(up_genes, expression$gene_name), ]
  up_exp <- as.data.frame(up_exp)

  duplicates <- up_exp$gene_name[which(duplicated(up_exp$gene_name))]

  if (!identical(duplicates, character(0))) {
    loc <- which(up_exp$gene_name == duplicates)
    strings <- paste0(up_exp$gene_name[loc], sep = "(", up_exp$gene_id[loc], ")")
    up_exp$gene_name[loc] <- strings
    up_genes <- up_exp$gene_name
  }

  up_exp <- up_exp[!is.na(up_exp$gene_name),]
  row.names(up_exp) <- up_exp$gene_name
  up_exp <- up_exp[,-c(1,2)]
  up_exp <- up_exp[,sample_groups[[comp]]]


  heatmap_up <- pheatmap(up_exp,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         fontsize_row = 40,
                         hm_colors(500),
                         scale= "row",
                         cellwidth = 150,
                         cellheight = 60,
                         angle_col=45,
                         fontsize_col=50)

  pdf(file.path(result_dir,paste0("UP_",groupName,"_heatmap.pdf")), width=50, height=top_DGEs)
  print(heatmap_up)
  dev.off()

  ##############################################################################

  down_genes <- df_down$gene_name
  down_genes <- down_genes[1:top_DGEs]


  down_exp <- expression[match(down_genes, expression$gene_name), ]

  down_exp <- as.data.frame(down_exp)

  duplicates <- down_exp$gene_name[which(duplicated(down_exp$gene_name))]

  if (!identical(duplicates, character(0))) {
    loc <- which(down_exp$gene_name == duplicates)
    strings <- paste0(down_exp$gene_name[loc], sep = "(", down_exp$gene_id[loc], ")")
    down_exp$gene_name[loc] <- strings
    down_genes <- down_exp$gene_name
  }

  down_exp <- down_exp[!is.na(down_exp$gene_name),]
  row.names(down_exp) <- down_exp$gene_name
  down_exp <- down_exp[,-c(1,2)]
  down_exp <- down_exp[,sample_groups[[comp]]]


  heatmap_down <- pheatmap(down_exp,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           fontsize_row = 40,
                           hm_colors(500),
                           scale= "row",
                           cellwidth = 150,
                           cellheight = 60,
                           angle_col=45,
                           fontsize_col=50)

  pdf(file.path(result_dir,paste0("DOWN_",groupName,"_heatmap.pdf")), width=50, height=top_DGEs)
  print(heatmap_down)
  dev.off()

  }

}
```

``` r
L <- list(c("S1", "S2", "S7", "S8"),
          c("S3", "S4", "S5", "S6"))

names(L) <- contrast_names

gsea.ora(deseq2_files, 
         proj_path = project_dir, 
         log2FC_cutoff = 1, 
         padj_cutoff = 0.05, 
         color_by = "p.adjust", 
         db = org.Mm.eg.db, 
         organism = "mmu", 
         minGSS = 3, 
         annots_df = annots_df,
         expression = expression,
         sample_groups = L,
         top_DGEs = 100)
```
