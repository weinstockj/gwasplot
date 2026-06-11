# Global variables to avoid NSE warnings
# utils::globalVariables(c(
#   "CHROM", "POS", "PVALUE", "gencode", "type", "gene_biotype", "chrom", 
#   "start", "end", "gene_id", "gene_name", "track", "y_center", "y_min", 
#   "y_max", "label_y", "distance"
# ))

#' Query ENCODE SCREEN cCREs for a genomic region
#' @param locus_chr Chromosome (e.g., "chr1")
#' @param locus_start Start position
#' @param locus_end End position  
#' @param assembly Genome assembly ("grch38" or "mm10")
#' @param biosample Optional biosample ID for cell-type-specific filtering (e.g., "GM12878_ENCDO000AAK")
#' @param cell_type Optional cell type name for filtering (e.g., "GM12878"). If provided, will be converted to biosample ID.
#' @return Data frame with cCRE information
#' @export
query_encode_ccres <- function(locus_chr, locus_start, locus_end, assembly = "grch38", biosample = NULL, cell_type = NULL) {  
  # Handle chromosome formatting - SCREEN expects "chr1" format
  if (!startsWith(locus_chr, "chr")) {
    locus_chr <- paste0("chr", locus_chr)
  }
  
  # Handle cell_type to biosample conversion if needed
  if (!is.null(cell_type) && is.null(biosample)) {
    # Get biosample ID from cell type name
    biosamples_df <- get_encode_biosamples(assembly)
    if (nrow(biosamples_df) > 0) {
      matching_biosample <- biosamples_df[grepl(cell_type, biosamples_df$display_name, ignore.case = TRUE), ]
      if (nrow(matching_biosample) > 0) {
        biosample <- matching_biosample$biosample[1]
        cli::cli_alert_info(paste0("Using biosample ID: ", biosample, " for cell type: ", cell_type))
      } else {
        warning("Cell type '", cell_type, "' not found. Available cell types can be viewed with get_encode_biosamples()")
        return(data.frame())
      }
    }
  }
    # Base GraphQL query
  if (is.null(biosample)) {
    # Query for general cCRE data
    query <- paste0('
    query {
      cCRESCREENSearch(
        assembly: "', assembly, '"
        coordinates: {
          chromosome: "', locus_chr, '"
          start: ', locus_start, '
          end: ', locus_end, '
        }
      ) {
        chrom
        start
        len
        pct
        ctcf_zscore
        dnase_zscore
        atac_zscore
        enhancer_zscore
        promoter_zscore
        info {
          accession
          isproximal
        }
      }
    }')
  } else {
    # Query for biosample-specific cCRE data
    query <- paste0('
    query {
      cCRESCREENSearch(
        assembly: "', assembly, '"
        coordinates: {
          chromosome: "', locus_chr, '"
          start: ', locus_start, '
          end: ', locus_end, '
        }
        cellType: "', biosample, '"
      ) {
        chrom
        start
        len
        pct
        ctspecific {
          ct
          ctcf_zscore
          dnase_zscore
          h3k4me3_zscore
          h3k27ac_zscore
          atac_zscore
        }
        info {
          accession
          isproximal
        }
      }    }')
  }
    # Make the GraphQL request
  tryCatch({
    # Prepare the request body exactly like PowerShell
    body <- list(query = query)
    body_json <- jsonlite::toJSON(body, auto_unbox = TRUE)
    
    response <- httr::POST(
      url = "https://screen.api.wenglab.org/graphql",
      body = body_json,
      httr::content_type("application/json"),
      encode = "raw"
    )
    
    if (httr::status_code(response) != 200) {
      cli::cli_warn("ENCODE SCREEN API request failed with status: ", httr::status_code(response))
      return(data.frame())
    }
    
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    result <- jsonlite::fromJSON(content)
    
    if (!is.null(result$errors)) {
      cli::cli_abort("ENCODE SCREEN API returned errors: ", paste(result$errors$message, collapse = "; "))
      return(data.frame())
    }    # Extract cCRE data
    ccres <- result$data$cCRESCREENSearch
    
    if (is.null(ccres) || length(ccres$chrom) == 0) {
      return(data.frame())
    }
      # Process the data - handle the structure properly
    ccres_df <- tibble::tibble(
      chrom = ccres$chrom,
      start = ccres$start,
      end = ccres$start + ccres$len,
      length = ccres$len,
      accession = ccres$info$accession,
      isproximal = ccres$info$isproximal,
      ccre_type = ccres$pct  # Use the pre-computed cCRE classification from the API
    )
    
    # Debug: Show cCRE type distribution
    if (nrow(ccres_df) > 0) {
      type_counts <- table(ccres_df$ccre_type)
      cli::cli_alert_info(paste0("cCRE types found: ", paste(names(type_counts), "=", type_counts, collapse = ", ")))
    }
    
    return(ccres_df)
    
  }, error = function(e) {
    cli::cli_abort("Error querying ENCODE SCREEN: ", e$message)
    return(data.frame())
  })
}

#' Assign genes to tracks to avoid overlaps
#' @param gene_df Data frame with gene information including start and end positions
#' @return Data frame with track assignments added
assign_gene_tracks <- function(gene_df) {
  if (nrow(gene_df) == 0) return(gene_df)
  
  # Sort genes by start position to process left to right
  gene_df <- gene_df[order(gene_df$start), ]
  n_genes <- nrow(gene_df)
  tracks <- rep(1, n_genes)
  
  if (n_genes == 1) {
    gene_df$track <- tracks
    return(gene_df)
  }
  
  # Buffer distance between genes (1kb)
  buffer <- 1000
  
  for (i in 2:n_genes) {
    # Find the lowest available track for gene i
    for (track in 1:i) {
      # Check if this gene overlaps with any gene already assigned to this track
      overlaps <- FALSE
      for (j in 1:(i-1)) {
        if (tracks[j] == track) {
          # Check for overlap with buffer
          if (gene_df$start[i] <= gene_df$end[j] + buffer) {
            overlaps <- TRUE
            break
          }
        }
      }
      if (!overlaps) {
        tracks[i] <- track
        break
      }
    }
  }
  
  gene_df$track <- tracks
  return(gene_df)
}

#' LocusZoom-style plot for GWAS results
#'
#' @param x A GWASFormatter object or a data.frame/tibble containing GWAS results.
#' @param locus_chr Chromosome of the locus to plot (e.g., 'chr1').
#' @param locus_start Start position of the locus (integer).
#' @param locus_end End position of the locus (integer).
#' @param include_ccres Logical, whether to include ENCODE SCREEN cCREs track. Default is FALSE.
#' @param ccre_biosample Optional biosample ID for biosample-specific cCRE filtering (e.g., "GM12878_ENCDO000AAK").
#' @param ccre_cell_type Optional cell type name for cCRE filtering (e.g., "GM12878"). Will be converted to biosample ID.
#' @param ... Additional arguments passed to methods.
#' @return A ggplot2 object with the locuszoom-style plot.
#' @export
locuszoom <- function(x, locus_chr, locus_start, locus_end, include_ccres = FALSE, ccre_biosample = NULL, ccre_cell_type = NULL, valid_biotype = "protein_coding", ...) {
  UseMethod("locuszoom")
}

#' @describeIn locuszoom Method for GWASFormatter objects
#' @export
locuszoom.GWASFormatter <- function(x, locus_chr, locus_start, locus_end, include_ccres = FALSE, ccre_biosample = NULL, ccre_cell_type = NULL, valid_biotype = "protein_coding", ...) {
  df <- x$data %>%
    dplyr::filter(CHROM == locus_chr, POS >= locus_start, POS <= locus_end) %>%
    dplyr::collect()
  locuszoom.data.frame(df, locus_chr, locus_start, locus_end, include_ccres, ccre_biosample, ccre_cell_type, valid_biotype = valid_biotype, ...)
}

#' @describeIn locuszoom Method for data.frame/tibble objects
#' @export
locuszoom.data.frame <- function(x, locus_chr, locus_start, locus_end, include_ccres = FALSE, ccre_biosample = NULL, ccre_cell_type = NULL, valid_biotype = "protein_coding", ...) {
  # Filter GWAS data for the locus
  df <- dplyr::filter(
    x, 
    CHROM == locus_chr, 
    POS >= locus_start, 
    POS <= locus_end
  )
  
  if (nrow(df) == 0) {
    stop("No variants found in the specified locus region")
  }
  # Get gene data for this region
  genes <- dplyr::filter(
    gencode,
    type == "gene",
    gene_biotype %in% .env[["valid_biotype"]],  # Use the valid_biotype parameter
    chrom == locus_chr,
    start <= locus_end,
    end >= locus_start
  )  # Get exon (CDS) data for this region


  exons <- dplyr::filter(
    gencode,
    type == "CDS",
    chrom == locus_chr,
    start <= locus_end,
    end >= locus_start
  )
    # Get cCRE data if requested
  ccres <- data.frame()
  if (include_ccres) {
    ccres <- query_encode_ccres(locus_chr, locus_start, locus_end, assembly = "grch38", biosample = ccre_biosample, cell_type = ccre_cell_type)
  }
  # Main association plot
  p_assoc <- ggplot2::ggplot(df, ggplot2::aes(x = POS, y = -log10(PVALUE))) +
    ggrastr::rasterize(
      ggplot2::geom_point(alpha = 0.7, color = "#2E4A87", size = 1.5),
      dpi = 300,
      dev = "ragg_png"
    ) +  # Dark blue points
    ggplot2::labs(
      title = paste0(locus_chr, ":", locus_start, "-", locus_end),
      x = "Position (Mb)",
      y = expression(-log[10](pvalue))
    ) +    
    ggplot2::scale_x_continuous(
      breaks = seq(min(df$POS), max(df$POS), length.out = 6),
      labels = scales::number(seq(min(df$POS), max(df$POS), length.out = 6) / 1e6, accuracy = 0.01),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))  # Add padding on both sides
    ) +    
    ggplot2::scale_y_continuous(
      breaks = {
        max_logp <- max(-log10(df$PVALUE))
        if (max_logp <= 30) {
          seq(0, ceiling(max_logp), by = 4)  # Every 4 for small ranges
        } else if (max_logp <= 50) {
          seq(0, ceiling(max_logp), by = 5)  # Every 10 for medium ranges
        } else {
          seq(0, ceiling(max_logp), by = 10)  # Every 20 for large ranges
        }
      },
      expand = ggplot2::expansion(mult = c(0.0, 0.1))
    ) +    cowplot::theme_cowplot(font_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      plot.margin = ggplot2::margin(t = 5, r = 10, b = 0, l = 5, unit = "pt")  # Reduced bottom margin from 10 to 5
    )
    # Create gene track plot
  if (nrow(genes) > 0) {
    # Assign genes to tracks to avoid overlaps
    genes <- assign_gene_tracks(genes)
    max_track <- max(genes$track)
    
    # Define heights: gene bodies are thin, exons are thick
    gene_body_height <- 0.15  # Thin gene body (introns)
    exon_height <- 0.4        # Thick exons
    track_spacing <- 0.6      # Space between tracks
    
    # Calculate y-positions for gene bodies (thin)
    genes <- genes %>%
      dplyr::mutate(
        y_center = (track - 1) * track_spacing,
        y_min = y_center - gene_body_height/2,
        y_max = y_center + gene_body_height/2,
        label_y = y_center + exon_height/2 + 0.1  # Position labels above the tallest part
      )
    
    # Update exon positions to be taller and centered on the same track
    if (nrow(exons) > 0) {
      exons <- exons %>%
        dplyr::left_join(
          genes %>% dplyr::select(gene_id, track, y_center),
          by = "gene_id"
        ) %>%
        dplyr::filter(!is.na(track)) %>%
        dplyr::mutate(
          y_min = y_center - exon_height/2,  # Taller exons
          y_max = y_center + exon_height/2
        )
    }
    # Create gene track plot
    p_genes <- ggplot2::ggplot() +
      # Plot gene bodies as thin rectangles (introns)
      ggplot2::geom_rect(
        data = genes,
        mapping = ggplot2::aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
        fill = "#cccccc", color = "#999999", alpha = 0.7
      )
    # Add exons as thicker rectangles on top
    if (nrow(exons) > 0) {
      p_genes <- p_genes +
        ggplot2::geom_rect(
          data = exons,
          mapping = ggplot2::aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
          fill = "#2c3e50", color = "#34495e", alpha = 0.9
        )
    }
    # Add gene labels positioned above genes in italics
    p_genes <- p_genes +
      ggplot2::geom_text(
        data = genes,
        mapping = ggplot2::aes(x = (start + end)/2, y = label_y, label = gene_name),
        size = 4.5,  # Increased from 3.9 to 4.5 for better readability
        hjust = 0.5,
        vjust = 0,
        fontface = "italic",
        color = "#2c3e50"
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-0.2, max_track * track_spacing + 0.1),  # Reduced bottom padding from -0.4 to -0.2
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::scale_x_continuous(
        limits = c(locus_start, locus_end),
        expand = c(0, 0)
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5, unit = "pt")
      )
  } else {
    # Create empty gene track if no genes found
    p_genes <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::scale_x_continuous(
        limits = c(locus_start, locus_end),
        expand = c(0, 0)
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5, unit = "pt")
      )
  }

  # Create cCRE track if requested
  if (include_ccres && nrow(ccres) > 0) {
    # Filter out inactive cCREs (InActive, noclass, TF) as they represent inactive/unclassified elements
    ccres_filtered <- ccres %>%
      dplyr::filter(ccre_type %in% c(
        "CA",
        "CA-H3K4me3",
        "CA-CTCF",
        "CA-TF",
        "PLS",
        "pELS",
        "dELS"
        # "TF"
      ))  # Keep only active cCRE types
    
    # Only proceed if we have active cCREs to plot
    if (nrow(ccres_filtered) > 0) {
      # Define colors for different cCRE types based on ENCODE SCREEN standards
      # Colors match the official ENCODE SCREEN visualization scheme
      ccre_colors <- c(
        "PLS" = "#FF0000",                 # Red - Promoter-like signatures
        "pELS" = "#FFA700",                # Orange - Proximal enhancer-like signatures
        "dELS" = "#FFA700",                # Orange - Distal enhancer-like signatures  
        "CA-H3K4me3" = "#FFCD00",          # Yellow - chromatin accessibility + H3K4me3
        "CA-CTCF" = "#00B0F0",             # Blue - chromatin accessibility + CTCF
        "CA" = "#00B050",                  # Green - chromatin accessibility only
        "CA-TF" = "#7030A0"                # Purple - chromatin accessibility + TF
      )
      
      p_ccres <- ggplot2::ggplot(ccres_filtered, ggplot2::aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0.2, fill = ccre_type)) +
      ggplot2::geom_rect(alpha = 0.8, linewidth = 0) +
      ggplot2::scale_fill_manual(values = ccre_colors, name = "cCRE Type") +
      ggplot2::scale_x_continuous(
        limits = c(locus_start, locus_end),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-0.4, 0.4),
        expand = ggplot2::expansion(mult = c(0, 0.01))
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 9),
        legend.key.size = ggplot2::unit(0.3, "cm"),
        legend.margin = ggplot2::margin(t = 2, b = 4, r = 10, l = 5),
        legend.box.spacing = ggplot2::unit(0.1, "cm")
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 0.3)))
    
      # Combine all three plots
      plot_heights <- c(5, 1, 1)  # Association, genes, cCREs - increased cCRE height
      cowplot::plot_grid(p_assoc, p_genes, p_ccres, ncol = 1, align = "v", rel_heights = plot_heights)
    } else {
      # No active cCREs to plot, combine just association and gene plots
      cowplot::plot_grid(p_assoc, p_genes, ncol = 1, align = "v", rel_heights = c(5, 1))
    }
  } else {
    # Combine just association and gene plots
    cowplot::plot_grid(p_assoc, p_genes, ncol = 1, align = "v", rel_heights = c(5, 1))
  }
}

#' Get available ENCODE SCREEN biosamples for a given assembly
#' @param assembly Genome assembly ("grch38" or "mm10")
#' @return Data frame with biosample metadata including cell type names
#' @export
get_encode_biosamples <- function(assembly = "grch38") {
  
  # GraphQL query to get biosample metadata
  query <- paste0('
  query {
    ccREBiosampleQuery(assembly: "', assembly, '") {
      biosamples {
        name
        ontology        
        lifeStage
        sampleType
        displayname
        dnase: experimentAccession(assay: "DNase")
        h3k4me3: experimentAccession(assay: "H3K4me3")
        h3k27ac: experimentAccession(assay: "H3K27ac")
        ctcf: experimentAccession(assay: "CTCF")
        atac: experimentAccession(assay: "ATAC")
      }
    }
  }')
  
  # Make the GraphQL request
  tryCatch({
    # Prepare the request body exactly like PowerShell
    body <- list(query = query)
    body_json <- jsonlite::toJSON(body, auto_unbox = TRUE)
    
    response <- httr::POST(
      url = "https://screen.api.wenglab.org/graphql",
      body = body_json,
      httr::content_type("application/json"),
      encode = "raw"
    )
    
    if (httr::status_code(response) != 200) {
      warning("ENCODE SCREEN API request failed with status: ", httr::status_code(response))
      return(data.frame())
    }
    
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    result <- jsonlite::fromJSON(content)
    
    if (!is.null(result$errors)) {
      warning("ENCODE SCREEN API returned errors: ", paste(result$errors$message, collapse = "; "))
      return(data.frame())
    }
    
    # Extract biosample data
    biosamples <- result$data$ccREBiosampleQuery$biosamples
    
    if (is.null(biosamples) || nrow(biosamples) == 0) {
      return(data.frame())
    }
    
    # Clean and organize the data
    biosamples_df <- tibble::tibble(
      biosample = biosamples$name,
      display_name = biosamples$displayname,
      ontology = biosamples$ontology,
      life_stage = biosamples$lifeStage,
      sample_type = biosamples$sampleType,
      has_dnase = !is.na(biosamples$dnase),
      has_h3k4me3 = !is.na(biosamples$h3k4me3),
      has_h3k27ac = !is.na(biosamples$h3k27ac),
      has_ctcf = !is.na(biosamples$ctcf),
      has_atac = !is.na(biosamples$atac),
      stringsAsFactors = FALSE
    )
    
    # Sort by display name for easier browsing
    biosamples_df <- biosamples_df[order(biosamples_df$display_name), ]
    
    return(biosamples_df)
    
  }, error = function(e) {
    warning("Error querying ENCODE SCREEN biosamples: ", e$message)
    return(data.frame())
  })
}

#' Search for available cell types and their biosample IDs
#' @param search_term Optional search term to filter cell types (case-insensitive)
#' @param assembly Genome assembly ("grch38" or "mm10")
#' @return Data frame with cell type information including biosample IDs
#' @export
search_cell_types <- function(search_term = NULL, assembly = "grch38") {
  biosamples_df <- get_encode_biosamples(assembly)
  
  if (nrow(biosamples_df) == 0) {
    return(data.frame())
  }
  
  # Filter by search term if provided
  if (!is.null(search_term)) {
    biosamples_df <- biosamples_df[grepl(search_term, biosamples_df$display_name, ignore.case = TRUE), ]
  }
  
  # Return relevant columns
  result <- biosamples_df[, c("biosample", "display_name", "ontology", "life_stage", "sample_type")]
  
  # Sort by display name
  result <- result[order(result$display_name), ]
  
  return(result)
}

#' Get biosample to accession mapping for a specific accession
#' @param accession cCRE accession number (e.g., "EH38E1234567")
#' @param assembly Genome assembly ("grch38" or "mm10") 
#' @return Data frame mapping biosample IDs to this accession
#' @export
get_biosample_for_accession <- function(accession, assembly = "grch38") {
  
  # GraphQL query to get biosample data for a specific accession
  query <- paste0('
  query {
    cCREQuery(
      assembly: "', assembly, '"
      accession: "', accession, '"
    ) {
      biosampleRankData {
        biosample
        rank
        rank_celltype
        zscore_h3k4me3
        zscore_h3k27ac
        zscore_ctcf
        zscore_dnase
        zscore_atac
      }
    }  }')
  
  # Make the GraphQL request
  tryCatch({
    # Prepare the request body exactly like PowerShell
    body <- list(query = query)
    body_json <- jsonlite::toJSON(body, auto_unbox = TRUE)
    
    response <- httr::POST(
      url = "https://screen.api.wenglab.org/graphql",
      body = body_json,
      httr::content_type("application/json"),
      encode = "raw"
    )
    
    if (httr::status_code(response) != 200) {
      warning("ENCODE SCREEN API request failed with status: ", httr::status_code(response))
      return(data.frame())
    }
    
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    result <- jsonlite::fromJSON(content)
    
    if (!is.null(result$errors)) {
      warning("ENCODE SCREEN API returned errors: ", paste(result$errors$message, collapse = "; "))
      return(data.frame())
    }
    
    # Extract biosample rank data
    rank_data <- result$data$cCREQuery$biosampleRankData
    
    if (is.null(rank_data) || nrow(rank_data) == 0) {
      return(data.frame())
    }
    
    # Convert to data frame and add metadata
    biosamples_df <- get_encode_biosamples(assembly)
    
    # Merge with biosample metadata
    result_df <- merge(
      rank_data, 
      biosamples_df[, c("biosample", "display_name", "ontology", "life_stage", "sample_type")],
      by = "biosample",
      all.x = TRUE
    )
    
    # Sort by rank (best ranks first)
    result_df <- result_df[order(result_df$rank), ]
    
    return(result_df)
    
  }, error = function(e) {
    warning("Error querying ENCODE SCREEN for accession: ", e$message)
    return(data.frame())
  })
}
