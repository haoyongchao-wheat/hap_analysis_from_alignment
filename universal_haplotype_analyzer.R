#!/usr/bin/env Rscript

# Universal Protein Haplotype Analyzer (R Version)
# 
# A comprehensive tool for analyzing protein haplotypes from amino acid alignment files.
# Generates high-quality visualizations and detailed reports suitable for publication.
# 
# Author: Haplotype Analysis Tool
# Version: 1.0.0
# License: MIT

# Load required libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(RColorBrewer)
  library(scales)
  library(gridExtra)
  library(svglite)
  library(optparse)
})

# Function to install missing packages
install_if_missing <- function(packages) {
  missing <- packages[!packages %in% installed.packages()[,"Package"]]
  if(length(missing) > 0) {
    cat("Installing missing packages:", paste(missing, collapse=", "), "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing, ask = FALSE)
  }
}

# Check and install required packages
required_packages <- c("Biostrings", "ggplot2", "dplyr", "readr", 
                      "jsonlite", "RColorBrewer", "scales", 
                      "gridExtra", "svglite", "optparse")
install_if_missing(required_packages)

# Universal Haplotype Analyzer Class
UniversalHaplotypeAnalyzer <- setRefClass(
  "UniversalHaplotypeAnalyzer",
  
  fields = list(
    fasta_file = "character",
    domain_file = "character",
    reference_id = "character",
    output_dir = "character",
    sequences = "list",
    reference_seq = "character",
    reference_name = "character",
    haplotypes = "list",
    polymorphisms = "data.frame",
    domains = "data.frame"
  ),
  
  methods = list(
    
    initialize = function(fasta_file, domain_file = NULL, reference_id = NULL, output_dir = ".") {
      "Initialize the analyzer"
      
      .self$fasta_file <- fasta_file
      .self$domain_file <- ifelse(is.null(domain_file), "", domain_file)
      .self$reference_id <- ifelse(is.null(reference_id), "", reference_id)
      .self$output_dir <- output_dir
      
      # Initialize data structures
      .self$sequences <- list()
      .self$reference_seq <- ""
      .self$reference_name <- ""
      .self$haplotypes <- list()
      .self$polymorphisms <- data.frame()
      .self$domains <- data.frame()
      
      # Create output directory
      if (!dir.exists(.self$output_dir)) {
        dir.create(.self$output_dir, recursive = TRUE)
      }
      
      # Validate inputs
      .self$validate_inputs()
    },
    
    validate_inputs = function() {
      "Validate input files"
      
      if (!file.exists(.self$fasta_file)) {
        stop(paste("FASTA file not found:", .self$fasta_file))
      }
      
      if (.self$domain_file != "" && !file.exists(.self$domain_file)) {
        stop(paste("Domain file not found:", .self$domain_file))
      }
    },
    
    load_sequences = function() {
      "Load sequences from FASTA file"
      
      cat("Loading sequences from", .self$fasta_file, "...\n")
      
      # Read FASTA file
      fasta_data <- readAAStringSet(.self$fasta_file)
      
      # Convert to character and remove stop codons and gaps
      sequence_order <- names(fasta_data)
      for (i in seq_along(fasta_data)) {
        seq_name <- names(fasta_data)[i]
        seq_str <- as.character(fasta_data[[i]])
        # Remove stop codons (*) and gaps (-)
        seq_str <- gsub("[*-]", "", seq_str)
        .self$sequences[[seq_name]] <- seq_str
      }
      
      if (length(.self$sequences) == 0) {
        stop("No sequences found in FASTA file")
      }
      
      # Set reference sequence
      if (.self$reference_id != "") {
        if (!(.self$reference_id %in% names(.self$sequences))) {
          stop(paste("Reference sequence not found:", .self$reference_id))
        }
        .self$reference_name <- .self$reference_id
      } else {
        # Use first sequence as reference
        .self$reference_name <- sequence_order[1]
      }
      
      .self$reference_seq <- .self$sequences[[.self$reference_name]]
      
      cat("Loaded", length(.self$sequences), "sequences\n")
      cat("Reference sequence:", .self$reference_name, 
          "(length:", nchar(.self$reference_seq), "aa)\n")
    },
    
    load_domains = function() {
      "Load domain information from file"
      
      if (.self$domain_file == "") {
        cat("No domain file provided, skipping domain analysis\n")
        return()
      }
      
      cat("Loading domain information from", .self$domain_file, "...\n")
      
      tryCatch({
        if (grepl("\\.json$", .self$domain_file, ignore.case = TRUE)) {
          .self$load_domains_json()
        } else {
          .self$load_domains_text()
        }
      }, error = function(e) {
        cat("Warning: Could not load domain file:", e$message, "\n")
        cat("Continuing without domain information\n")
        .self$domains <- data.frame()
      })
    },
    
    load_domains_json = function() {
      "Load domains from JSON file"
      
      domain_data <- fromJSON(.self$domain_file)
      
      domain_list <- list()
      for (domain_name in names(domain_data)) {
        info <- domain_data[[domain_name]]
        domain_list[[length(domain_list) + 1]] <- list(
          name = domain_name,
          start = info$start,
          end = info$end,
          color = ifelse(is.null(info$color), NA, info$color)
        )
      }
      
      .self$domains <- do.call(rbind, lapply(domain_list, data.frame, stringsAsFactors = FALSE))
    },
    
    load_domains_text = function() {
      "Load domains from text file"
      
      lines <- readLines(.self$domain_file)
      domain_list <- list()
      
      for (i in seq_along(lines)) {
        line <- trimws(lines[i])
        if (line == "" || startsWith(line, "#")) next
        
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) < 3) {
          cat("Warning: Invalid line", i, "in domain file:", line, "\n")
          next
        }
        
        tryCatch({
          domain_list[[length(domain_list) + 1]] <- list(
            name = parts[1],
            start = as.numeric(parts[2]),
            end = as.numeric(parts[3]),
            color = ifelse(length(parts) > 3, parts[4], NA)
          )
        }, error = function(e) {
          cat("Warning: Invalid coordinates in line", i, ":", line, "\n")
        })
      }
      
      if (length(domain_list) > 0) {
        .self$domains <- do.call(rbind, lapply(domain_list, data.frame, stringsAsFactors = FALSE))
      } else {
        .self$domains <- data.frame()
      }
    },
    
    generate_colors = function(n_colors) {
      "Generate n distinct colors"
      
      if (n_colors <= 8) {
        colors <- brewer.pal(max(3, n_colors), "Set2")[1:n_colors]
      } else {
        # Generate colors using HSV space
        hues <- seq(0, 1, length.out = n_colors + 1)[1:n_colors]
        colors <- hsv(hues, s = 0.7, v = 0.8)
      }
      
      return(colors)
    },
    
    assign_domain_colors = function() {
      "Assign colors to domains"
      
      if (nrow(.self$domains) == 0) return()
      
      # Check which domains need colors
      domains_without_colors <- which(is.na(.self$domains$color))
      
      if (length(domains_without_colors) > 0) {
        new_colors <- .self$generate_colors(length(domains_without_colors))
        .self$domains$color[domains_without_colors] <- new_colors
      }
      
      cat("Assigned colors to", nrow(.self$domains), "domains\n")
    },
    
    find_polymorphisms = function() {
      "Find polymorphic sites relative to reference sequence"
      
      cat("Identifying polymorphic sites...\n")
      
      ref_len <- nchar(.self$reference_seq)
      poly_list <- list()
      
      for (pos in 1:ref_len) {
        ref_aa <- substr(.self$reference_seq, pos, pos)
        variants <- c()
        
        for (seq_name in names(.self$sequences)) {
          if (seq_name == .self$reference_name) next
          
          sequence <- .self$sequences[[seq_name]]
          if (pos <= nchar(sequence)) {
            seq_aa <- substr(sequence, pos, pos)
            if (seq_aa != ref_aa) {
              variants <- c(variants, seq_aa)
            }
          }
        }
        
        if (length(variants) > 0) {
          poly_list[[length(poly_list) + 1]] <- list(
            position = pos,
            reference = ref_aa,
            variants = paste(unique(variants), collapse = ",")
          )
        }
      }
      
      if (length(poly_list) > 0) {
        .self$polymorphisms <- do.call(rbind, lapply(poly_list, data.frame, stringsAsFactors = FALSE))
      } else {
        .self$polymorphisms <- data.frame(position = integer(), 
                                         reference = character(), 
                                         variants = character())
      }
      
      cat("Found", nrow(.self$polymorphisms), "polymorphic sites\n")
    },
    
    classify_haplotypes = function() {
      "Classify and merge haplotypes"
      
      cat("Classifying haplotypes...\n")
      
      haplotype_groups <- list()
      
      for (seq_name in names(.self$sequences)) {
        if (seq_name == .self$reference_name) next
        
        sequence <- .self$sequences[[seq_name]]
        
        # Generate haplotype signature
        signature_parts <- c()
        for (i in 1:nrow(.self$polymorphisms)) {
          pos <- .self$polymorphisms$position[i]
          if (pos <= nchar(sequence)) {
            aa <- substr(sequence, pos, pos)
            signature_parts <- c(signature_parts, paste0(pos, ":", aa))
          } else {
            signature_parts <- c(signature_parts, paste0(pos, ":-"))
          }
        }
        
        signature <- paste(signature_parts, collapse = "|")
        
        if (signature %in% names(haplotype_groups)) {
          haplotype_groups[[signature]] <- c(haplotype_groups[[signature]], seq_name)
        } else {
          haplotype_groups[[signature]] <- seq_name
        }
      }
      
      # Assign names to haplotype groups
      haplotype_counter <- 1
      for (signature in names(haplotype_groups)) {
        seq_ids <- haplotype_groups[[signature]]
        
        if (length(seq_ids) == 1) {
          haplotype_name <- seq_ids[1]
        } else {
          haplotype_name <- paste0("Haplotype_", haplotype_counter)
          haplotype_counter <- haplotype_counter + 1
        }
        
        representative_seq <- .self$sequences[[seq_ids[1]]]
        .self$haplotypes[[haplotype_name]] <- list(
          sequence = representative_seq,
          members = seq_ids,
          signature = signature
        )
      }
      
      cat("Identified", length(.self$haplotypes), "distinct haplotypes\n")
      for (name in names(.self$haplotypes)) {
        info <- .self$haplotypes[[name]]
        cat("  ", name, ":", length(info$members), "member(s)\n")
      }
    },
    
    generate_report = function() {
      "Generate comprehensive analysis report"
      
      report_file <- file.path(.self$output_dir, "haplotype_analysis_report.txt")
      
      cat("Generating analysis report:", report_file, "\n")
      
      sink(report_file)
      
      cat("PROTEIN HAPLOTYPE ANALYSIS REPORT\n")
      cat(paste(rep("=", 50), collapse = ""), "\n\n")
      
      # Basic information
      cat("BASIC INFORMATION\n")
      cat(paste(rep("-", 20), collapse = ""), "\n")
      cat("Input file:", basename(.self$fasta_file), "\n")
      cat("Reference sequence:", .self$reference_name, "\n")
      cat("Reference length:", nchar(.self$reference_seq), "amino acids\n")
      cat("Total sequences analyzed:", length(.self$sequences), "\n")
      cat("Domain file:", ifelse(.self$domain_file == "", "None", basename(.self$domain_file)), "\n")
      cat("Number of domains:", nrow(.self$domains), "\n\n")
      
      # Polymorphism summary
      cat("POLYMORPHISM ANALYSIS\n")
      cat(paste(rep("-", 25), collapse = ""), "\n")
      cat("Total polymorphic sites:", nrow(.self$polymorphisms), "\n")
      if (nrow(.self$polymorphisms) > 0) {
        density <- nrow(.self$polymorphisms) / nchar(.self$reference_seq) * 100
        cat("Polymorphism density:", sprintf("%.2f%%", density), "\n\n")
        
        cat("Polymorphic sites details:\n")
        for (i in 1:nrow(.self$polymorphisms)) {
          poly <- .self$polymorphisms[i, ]
          cat(sprintf("  %2d. Position %4d: %s -> %s\n", 
                     i, poly$position, poly$reference, poly$variants))
        }
        cat("\n")
      } else {
        cat("\n")
      }
      
      # Haplotype summary
      cat("HAPLOTYPE ANALYSIS\n")
      cat(paste(rep("-", 20), collapse = ""), "\n")
      cat("Number of distinct haplotypes:", length(.self$haplotypes), "\n\n")
      
      for (i in seq_along(.self$haplotypes)) {
        name <- names(.self$haplotypes)[i]
        info <- .self$haplotypes[[name]]
        
        cat("Haplotype", i, ":", name, "\n")
        cat("  Sequence length:", nchar(info$sequence), "amino acids\n")
        cat("  Number of members:", length(info$members), "\n")
        cat("  Members:", paste(info$members, collapse = ", "), "\n")
        
        # Calculate differences from reference
        differences <- 0
        for (j in 1:nrow(.self$polymorphisms)) {
          pos <- .self$polymorphisms$position[j]
          if (pos <= nchar(info$sequence) && pos <= nchar(.self$reference_seq)) {
            if (substr(info$sequence, pos, pos) != substr(.self$reference_seq, pos, pos)) {
              differences <- differences + 1
            }
          }
        }
        
        cat("  Differences from reference:", differences, "sites\n")
        
        # Check for truncation
        if (nchar(info$sequence) < nchar(.self$reference_seq)) {
          cat("  Status: Truncated at position", nchar(info$sequence), "\n")
        } else {
          cat("  Status: Full length\n")
        }
        cat("\n")
      }
      
      # Domain information
      if (nrow(.self$domains) > 0) {
        cat("DOMAIN INFORMATION\n")
        cat(paste(rep("-", 18), collapse = ""), "\n")
        for (i in 1:nrow(.self$domains)) {
          domain <- .self$domains[i, ]
          cat("  ", domain$name, ":", domain$start, "-", domain$end, 
              "aa (color:", domain$color, ")\n")
        }
        cat("\n")
      }
      
      # Sequence statistics
      cat("SEQUENCE STATISTICS\n")
      cat(paste(rep("-", 19), collapse = ""), "\n")
      lengths <- sapply(.self$sequences, nchar)
      cat("  Mean length:", sprintf("%.1f", mean(lengths)), "aa\n")
      cat("  Length range:", min(lengths), "-", max(lengths), "aa\n")
      cat("  Standard deviation:", sprintf("%.1f", sd(lengths)), "aa\n")
      
      # Length distribution
      length_counts <- table(lengths)
      cat("\n  Length distribution:\n")
      for (length_val in names(length_counts)) {
        count <- length_counts[length_val]
        cat("    ", length_val, "aa:", count, "sequence(s)\n")
      }
      
      sink()
      
      return(report_file)
    },
    
    create_png_visualization = function() {
      "Create PNG visualization"
      
      png_file <- file.path(.self$output_dir, "haplotype_diagram.png")
      
      cat("Creating PNG visualization:", png_file, "\n")
      
      # Prepare data for plotting
      plot_data <- .self$prepare_plot_data()
      
      # Create the plot
      p <- .self$create_ggplot(plot_data)
      
      # Save as PNG
      ggsave(png_file, plot = p, width = 16, height = max(8, length(.self$haplotypes) * 0.8), 
             dpi = 300, units = "in")
      
      return(png_file)
    },
    
    create_svg_visualization = function() {
      "Create editable SVG visualization"
      
      svg_file <- file.path(.self$output_dir, "haplotype_diagram_editable.svg")
      
      cat("Creating editable SVG visualization:", svg_file, "\n")
      
      # Prepare data for plotting
      plot_data <- .self$prepare_plot_data()
      
      # Create the plot
      p <- .self$create_ggplot(plot_data)
      
      # Save as SVG
      ggsave(svg_file, plot = p, width = 16, height = max(8, length(.self$haplotypes) * 0.8), 
             units = "in")
      
      return(svg_file)
    },
    
    prepare_plot_data = function() {
      "Prepare data for plotting"
      
      seq_length <- nchar(.self$reference_seq)
      
      # Sequence data
      seq_data <- data.frame(
        name = character(),
        start = numeric(),
        end = numeric(),
        y_pos = numeric(),
        type = character(),
        stringsAsFactors = FALSE
      )
      
      # Reference sequence
      seq_data <- rbind(seq_data, data.frame(
        name = .self$reference_name,
        start = 1,
        end = seq_length,
        y_pos = length(.self$haplotypes) + 1,
        type = "reference",
        stringsAsFactors = FALSE
      ))
      
      # Haplotype sequences
      y_pos <- length(.self$haplotypes)
      for (name in names(.self$haplotypes)) {
        info <- .self$haplotypes[[name]]
        display_name <- name
        if (length(info$members) > 1) {
          display_name <- paste0(name, " (", length(info$members), " members)")
        }
        
        seq_data <- rbind(seq_data, data.frame(
          name = display_name,
          start = 1,
          end = nchar(info$sequence),
          y_pos = y_pos,
          type = "haplotype",
          stringsAsFactors = FALSE
        ))
        y_pos <- y_pos - 1
      }
      
      # Domain data
      domain_data <- data.frame(
        name = character(),
        domain_name = character(),
        start = numeric(),
        end = numeric(),
        y_pos = numeric(),
        color = character(),
        stringsAsFactors = FALSE
      )
      
      if (nrow(.self$domains) > 0) {
        # Add domains for reference
        for (i in 1:nrow(.self$domains)) {
          domain <- .self$domains[i, ]
          if (domain$start <= seq_length) {
            domain_end <- min(domain$end, seq_length)
            domain_data <- rbind(domain_data, data.frame(
              name = .self$reference_name,
              domain_name = domain$name,
              start = domain$start,
              end = domain_end,
              y_pos = length(.self$haplotypes) + 1,
              color = domain$color,
              stringsAsFactors = FALSE
            ))
          }
        }
        
        # Add domains for haplotypes
        y_pos <- length(.self$haplotypes)
        for (name in names(.self$haplotypes)) {
          info <- .self$haplotypes[[name]]
          actual_length <- nchar(info$sequence)
          
          for (i in 1:nrow(.self$domains)) {
            domain <- .self$domains[i, ]
            if (domain$start <= actual_length) {
              domain_end <- min(domain$end, actual_length)
              if (domain_end > domain$start) {
                domain_data <- rbind(domain_data, data.frame(
                  name = name,
                  domain_name = domain$name,
                  start = domain$start,
                  end = domain_end,
                  y_pos = y_pos,
                  color = domain$color,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
          y_pos <- y_pos - 1
        }
      }
      
      # Polymorphism data
      poly_data <- data.frame(
        position = numeric(),
        y_pos = numeric(),
        stringsAsFactors = FALSE
      )
      
      if (nrow(.self$polymorphisms) > 0) {
        y_pos <- length(.self$haplotypes)
        for (name in names(.self$haplotypes)) {
          info <- .self$haplotypes[[name]]
          actual_length <- nchar(info$sequence)
          
          for (i in 1:nrow(.self$polymorphisms)) {
            pos <- .self$polymorphisms$position[i]
            if (pos <= actual_length && pos <= nchar(.self$reference_seq)) {
              if (substr(info$sequence, pos, pos) != substr(.self$reference_seq, pos, pos)) {
                poly_data <- rbind(poly_data, data.frame(
                  position = pos,
                  y_pos = y_pos,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
          y_pos <- y_pos - 1
        }
      }
      
      return(list(
        sequences = seq_data,
        domains = domain_data,
        polymorphisms = poly_data,
        seq_length = seq_length
      ))
    },
    
    create_ggplot = function(plot_data) {
      "Create ggplot visualization"
      
      seq_length <- plot_data$seq_length
      
      # Base plot
      p <- ggplot() +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 1),
          plot.margin = margin(20, 20, 20, 100)
        ) +
        scale_x_continuous(
          name = "Amino acid position",
          limits = c(0, seq_length),
          breaks = seq(0, seq_length, by = max(100, seq_length %/% 10)),
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          limits = c(0.5, length(.self$haplotypes) + 1.5),
          expand = c(0, 0)
        )
      
      # Add sequence bars
      p <- p + geom_rect(
        data = plot_data$sequences,
        aes(xmin = start, xmax = end, ymin = y_pos - 0.3, ymax = y_pos + 0.3),
        fill = ifelse(plot_data$sequences$type == "reference", "#ECF0F1", "white"),
        color = "#2C3E50",
        size = 1
      )
      
      # Add domains
      if (nrow(plot_data$domains) > 0) {
        p <- p + geom_rect(
          data = plot_data$domains,
          aes(xmin = start, xmax = end, ymin = y_pos - 0.3, ymax = y_pos + 0.3, fill = color),
          alpha = 0.8,
          color = "white",
          size = 0.5
        ) +
        scale_fill_identity()
      }
      
      # Add polymorphisms
      if (nrow(plot_data$polymorphisms) > 0) {
        p <- p + geom_segment(
          data = plot_data$polymorphisms,
          aes(x = position, xend = position, y = y_pos - 0.4, yend = y_pos + 0.4),
          color = "#E74C3C",
          size = 1.5
        )
      }
      
      # Add sequence labels
      label_data <- plot_data$sequences
      p <- p + geom_text(
        data = label_data,
        aes(x = -seq_length * 0.05, y = y_pos, label = name),
        hjust = 1,
        vjust = 0.5,
        size = 4,
        fontface = "bold"
      )
      
      # Add truncation lines
      truncated_seqs <- plot_data$sequences[plot_data$sequences$end < seq_length, ]
      if (nrow(truncated_seqs) > 0) {
        p <- p + geom_segment(
          data = truncated_seqs,
          aes(x = end, xend = seq_length, y = y_pos, yend = y_pos),
          linetype = "dashed",
          color = "#95A5A6",
          size = 1
        ) +
        geom_segment(
          data = truncated_seqs,
          aes(x = end, xend = end, y = y_pos - 0.2, yend = y_pos + 0.2),
          color = "#E74C3C",
          size = 2
        )
      }
      
      return(p)
    },
    
    run_analysis = function() {
      "Run complete haplotype analysis"
      
      cat("Starting Universal Protein Haplotype Analysis...\n")
      cat(paste(rep("=", 50), collapse = ""), "\n")
      
      tryCatch({
        # Load data
        .self$load_sequences()
        .self$load_domains()
        .self$assign_domain_colors()
        
        # Analyze
        .self$find_polymorphisms()
        .self$classify_haplotypes()
        
        # Generate outputs
        report_file <- .self$generate_report()
        png_file <- .self$create_png_visualization()
        svg_file <- .self$create_svg_visualization()
        
        cat("\n", paste(rep("=", 50), collapse = ""), "\n")
        cat("Analysis completed successfully!\n")
        cat("\nGenerated files:\n")
        cat("  - Report:", report_file, "\n")
        cat("  - PNG visualization:", png_file, "\n")
        cat("  - Editable SVG:", svg_file, "\n")
        
        return(list(
          report = report_file,
          png = png_file,
          svg = svg_file
        ))
        
      }, error = function(e) {
        cat("\nError during analysis:", e$message, "\n")
        stop(e)
      })
    }
  )
)

# Main function
main <- function() {
  # Command line options
  option_list <- list(
    make_option(c("-d", "--domains"), type = "character", default = NULL,
                help = "Domain information file (tab-separated or JSON)"),
    make_option(c("-r", "--reference"), type = "character", default = NULL,
                help = "Reference sequence ID (default: first sequence)"),
    make_option(c("-o", "--output"), type = "character", default = ".",
                help = "Output directory (default: current directory)"),
    make_option(c("-h", "--help"), action = "store_true", default = FALSE,
                help = "Show this help message")
  )
  
  parser <- OptionParser(
    usage = "Usage: %prog [options] fasta_file",
    option_list = option_list,
    description = "Universal Protein Haplotype Analyzer (R Version)",
    epilogue = paste(
      "Examples:",
      "  # Basic analysis",
      "  Rscript universal_haplotype_analyzer.R sequences.fasta",
      "",
      "  # With domain information",
      "  Rscript universal_haplotype_analyzer.R -d domains.txt sequences.fasta",
      "",
      "  # Specify reference sequence and output directory",
      "  Rscript universal_haplotype_analyzer.R -r REF_SEQ -o results/ sequences.fasta",
      sep = "\n"
    )
  )
  
  # Parse arguments
  args <- parse_args(parser, positional_arguments = 1)
  opt <- args$options
  
  if (length(args$args) != 1) {
    print_help(parser)
    quit(status = 1)
  }
  
  fasta_file <- args$args[1]
  
  tryCatch({
    analyzer <- UniversalHaplotypeAnalyzer$new(
      fasta_file = fasta_file,
      domain_file = opt$domains,
      reference_id = opt$reference,
      output_dir = opt$output
    )
    
    analyzer$run_analysis()
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}