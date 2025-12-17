library(DT)
library(shinythemes)
library(shiny)
library(shinydashboard)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(biomaRt)
library(Biostrings)
library(tidyverse)
library(cluster)
library(rlang)
library(dplyr)
library(purrr)
library(shinyBS)
library(openxlsx)

source("../tools/GGGenome_functions.R")
source("../tools/RNaseH_script.R")
source("../tools/offtarget.R")

# ----------------------------------- working directory ------------------------
working_dir = getwd()
setwd(working_dir)

print("work directory set")

function(input, output, session) {
  
  # ----------------------------------- Notifications UI -----------------------
  # Notification stays until clicked away
  showNotification(
    "Finished loading",
    type = "default",
    duration = NULL,
    closeButton = TRUE
  )

  t1 <- Sys.time()
  
  observeEvent(input$run_button, {
    showNotification(
      "Script started",
      type = "default",
      duration = NULL,
      closeButton = TRUE
    )
  
  # ----------------------------------- Functions ----------------------------
  # Make the tox score function
  
  calculate_acute_neurotox <- function(xx) {
    # make sure input is in character format
    xx <- as.character(xx)
    
    # count the number of each nucleotide
    lf <- function(x) {
      x <- tolower(x)
      x <- strsplit(x, "")[[1]]
      x <- table(factor(x, levels = c("a", "c", "t", "g")))
      return(x)
    }
    cnt_nt <- as.data.frame(t(sapply(xx, lf)))
    # count number of nucleotides from the 3'-end untill it finds the first g
    gfree3 <- function(x) {
      x <- tolower(x)
      x <- strsplit(x, "")[[1]]
      tfg <- x == "g"
      if (sum(tfg) == 0) {
        l3 <- NA
      } else {
        posg <- c(1:length(x))[tfg]
        l3 <- length(x) - max(posg)
      }
      return(l3)
    }
    cnt_gfree3 <- sapply(xx, gfree3)
    cnt_gfree3[cnt_gfree3 > 20] <- 20      #Set max to 20
    cnt_gfree3[is.na(cnt_gfree3)] <- 20  #Set no g in ASO to 20
    ## Calculate final score based on trained parameters and return result
    calc_out <- round(
      136.0430 - 3.1263 * cnt_nt$a - 5.1100 * cnt_nt$c -
        4.7217 * cnt_nt$t - 10.1264 * cnt_nt$g + 1.3577 *
        cnt_gfree3,
      1
    )
    
    return(as.numeric(calc_out))
  }
  
    
  if (input$linux_input == TRUE) {
    # 2.5.1 Predict accessibility for an RNA target
    RNAplfold_R = function(seq.char,
                           L.in = 40,
                           W.in = 80,
                           u.in = 16) {
      cmmnd2 = paste("RNAplfold -L", L.in, "-W", W.in, "-u", u.in)
      seq.char = as.character(seq.char)
      cat(seq.char, file = paste("|", cmmnd2, sep = ""))
      acc.tx = read.delim(
        "plfold_lunp",
        as.is = T,
        skip = 2,
        header = F,
        row.names = 1
      )
      acc.tx = acc.tx[, colSums(is.na(acc.tx)) != nrow(acc.tx)]
      colnames(acc.tx) = 1:ncol(acc.tx)
      file.remove("plfold_lunp")
      file.remove("plfold_dp.ps")
      return(acc.tx)
    }
    
    # 2.5.3 Predict duplex formation and self folding of oligonucleotides
    RNAduplex_R = function(seqs) {
      sys_cmd = system('RNAduplex',
                       input = c(seqs, seqs),
                       intern = TRUE)
      as.numeric(regmatches(sys_cmd, regexpr("-?\\d+\\.\\d+", sys_cmd)))
    }


    # --- Jims data location with docker --- 
    # txdb_hsa <- loadDb("/opt/ASOstool-v2/txdb_hsa_biomart.db")
    
    # --- Harrys data location with script ---
    txdb_hsa <- loadDb("../Data/txdb_hsa.db")
    
    # ----------------------------------- milestone 1 --------------------------
    print("milestone1")

    
    RNAselffold_R = function (seqs) {
      output = system('RNAfold --noPS',
                      input = c(seqs),
                      intern = T)
      output = unlist(strsplit(output[grepl('[0-9]', output)], '[(]'))
      as.double(gsub(' |[)]', '', output[grepl('[0-9]', output)]))
    }
  }
  
  filter_function <- function(df,
                              valueX,
                              columname,
                              operator_string) {
    switch(
      operator_string,
      "==" = {
        # Code for when my_string is "=="
        filtered_df <- df %>% filter(.data[[columname]] == valueX)
        return(filtered_df)
      },
      "!=" = {
        # Code for when my_string is "!="
        filtered_df <- df %>% filter(.data[[columname]] != valueX)
        return(filtered_df)
      },
      "<" = {
        # Code for when my_string is "<"
        filtered_df <- df %>% filter(.data[[columname]] < valueX)
        return(filtered_df)
      },
      ">" = {
        # Code for when my_string is ">"
        filtered_df <- df %>% filter(.data[[columname]] > valueX)
        return(filtered_df)
      },
      "<=" = {
        # Code for when my_string is "<="
        filtered_df <- df %>% filter(.data[[columname]] <= valueX)
        return(filtered_df)
      },
      ">=" = {
        # Code for when my_string is ">="
        filtered_df <- df %>% filter(.data[[columname]] >= valueX)
        return(filtered_df)
      },
      {
        # Default case
        return(df)
      }
    )
  }
  # ----------------------------------- Data setup ---------------------------
  # Store all human pre-mRNA sequences
  
  # path = getwd()
  
  # txdb_hsa <- loadDb("txdb_hsa_biomart.db")
  getwd()
  # --- Harrys data location with script ---
  txdb_hsa <- loadDb("../txdb_hsa_biomart.db")
  txdb_hsa <- loadDb("../Data/txdb_hsa.db")
    
  # ----------------------------------- milestone 1 --------------------------
  print("milestone1: loaded human database object")
  
  # Extract the genes
  gdb_hsa <- genes(txdb_hsa)
  
  # Define the chromosomes to keep
  chr_to_keep <- c(as.character(1:22), 'X', 'Y', 'MT')
  
  # Filtert Hsapiens zodat het alleen de aangegeven chromosomen pakt
  Hsapiens@user_seqnames <- setNames(chr_to_keep, chr_to_keep)
  Hsapiens@seqinfo <- Hsapiens@seqinfo[chr_to_keep]
  
  # Subset genes to keep only those on specified chromosomes
  gdb_hsa <- gdb_hsa[seqnames(gdb_hsa) %in% chr_to_keep]
  
  # ----------------------------------- milestone 2 --------------------------
  print("milestone2: Subsetted genes from specified chromosomes")
  
  # Get the sequences *
  HS <- getSeq(Hsapiens, gdb_hsa)
    
  # ----------------------------------- milestone 3 --------------------------
  print("milestone3: Saved human gene sequences")
  # Target collect the pre-mRNA sequence
  
  # Define wanted Ensembl ID
  ensembl_ID = input$ensemble_id_input
  
  # Retrieve a specific RNA target using the Ensembl ID
  RNA_target = HS[names(HS) == ensembl_ID]
  
  
  # ----------------------------------- milestone 4 --------------------------
  print("milestone4: Saved input RNA target (ENSEMBL) information: \n")
  print(RNA_target)
  #filters on ensembl ID
  target_ranges = gdb_hsa[names(gdb_hsa) == ensembl_ID]

  
  # Extracts the chromosome name for target range,extracts the start and end
  # Position of the genomic region and note from which strand it is.
  # Positive strand 1 ('+'), negative strand -1 ('-'), or unspecified 0 ('*').
  chr_coord = c(
    chr = as.numeric(as.character(seqnames(target_ranges))),
    start = start(target_ranges),
    end = end(target_ranges),
    strand = ifelse(strand(target_ranges) == "+", 1, -1)
  )
  
  # ----------------------------------- milestone 5 --------------------------
  print("milestone5: Extracted chromosome coordinates of RNA target")
    
  # Set the width of rna_target as value L
  l = width(RNA_target)
  
  # Note length of subsequence
  oligo_lengths = input$oligo_length_range[1]:input$oligo_length_range[2]
  
  # Construct the dataframe
  target_annotation = lapply(oligo_lengths, function(i) {
    tibble(start = 1:(l - i + 1), length = i)
  }) %>%
    bind_rows() %>%
    mutate(end = start + length - 1)
  
  # Adds "name" sequence to the dataframe
  target_regions = DNAStringSet(RNA_target[[1]],
                                start = target_annotation$start,
                                width = target_annotation$length)
  names(target_regions) = as.character(target_regions)
  
  target_annotation$name = names(target_regions)
  
  # ----------------------------------- milestone 6 --------------------------
  print("milestone6: Enumerated all possible ASO target sequences")
  prefilter <- nrow(target_annotation)
  
  target_annotation <- target_annotation %>%
    filter(!grepl("^C", name))
  
  postfilter <- nrow(target_annotation)
  removed <- prefilter - postfilter
  
  print("Filtering Oligo sequences ending with G.")
  print(paste0("Rows before filtering: ", prefilter))
  print(paste0("Rows after filtering: ", postfilter))
  print(paste0("Filtering removed ", removed, " possible ASOs."))

  # ----------------------------------- milestone 7 --------------------------
  print("milestone 7: Prefiltered Oligo sequences ending with G")
  
  if (input$linux_input == TRUE) {
    # 3.4 Estimate Transcript Accessibility for the RNA Target at Single-Nucleotide Resolution
    accessibility = RNAplfold_R(RNA_target, u.in = max(oligo_lengths)) %>%
      as_tibble() %>%
      mutate(end = 1:l) %>%
      gather(length, accessibility, -end) %>%
      mutate(length = as.double(length))
    
    target_annotation = left_join(target_annotation, accessibility, by =
                                    c('length', 'end'))
    print("filtering accesibility")
    
    prefilter <- nrow(target_annotation)
    accessibility_cutoff <- input$numeric_input_c
    target_annotation <- target_annotation %>%
      filter(accessibility > accessibility_cutoff)
    
    postfilter <- nrow(target_annotation)
    removed <- prefilter - postfilter
    
    print(paste0("Filtering Oligo sequences with accissibility score > ", accessibility_cutoff))
    print(paste0("Rows before filtering: ", prefilter))
    print(paste0("Rows after filtering: ", postfilter))
    print(paste0("Filtering removed ", removed, " possible ASOs."))
    }
  
  # ----------------------------------- milestone 8 --------------------------
  print("milestone 8: Calculated ViennaRNA accessibility score and filtering")
  
  nucleobase_seq = reverseComplement(target_regions)
  
  # Voeg ze toe aan uni_tar gekoppeld aan name
  target_annotation$oligo_seq = as.character(nucleobase_seq[target_annotation$name])
  
  # Toxicity score acute neurotox score (desired >70)
  target_annotation$tox_score = calculate_acute_neurotox(target_annotation$oligo_seq)
  
  prefilter <- nrow(target_annotation)
  tox_score_cutoff <- input$numeric_input_e
  target_annotation <- target_annotation %>%
    filter(tox_score > tox_score_cutoff)
  
  postfilter <- nrow(target_annotation)
  removed <- prefilter - postfilter
  
  print(paste0("Filtering Oligo sequences with toxicity score > ", tox_score_cutoff))
  print(paste0("Rows before filtering: ", prefilter))
  print(paste0("Rows after filtering: ", postfilter))
  print(paste0("Filtering removed ", removed, " possible ASOs."))
  
  # ----------------------------------- milestone 9 --------------------------
  print("milestone 9: Calculated toxicity score and filtering")
  
  # Bereken aantal "cg"
  target_annotation$CGs = (target_annotation$length -
                   nchar(gsub('CG', '', target_annotation$name))) / 2
  
  # ----------------------------------- milestone 10 --------------------------
  print("milestone 10: Calculated CG motifs")
  
  
  # Define the marts for mmusculus and hsapiens
  martHS = useEnsembl(biomart="ensembl",
                      dataset="hsapiens_gene_ensembl")
  if (input$Conserved_input == TRUE) {
  martMM = useEnsembl(biomart="ensembl",
                      dataset="mmusculus_gene_ensembl")
  
  # Get the orthologous Ensembl gene for the provided human Ensembl ID
  ortho_ENS = getBM(attributes = "mmusculus_homolog_ensembl_gene",
                    filters = "ensembl_gene_id",
                    values = ensembl_ID, mart = martHS,
                    bmHeader = FALSE)
  
  # ----------------------------------- milestone 11 --------------------------
  print("milestone 11: Get mouse ortholog data for genes")
  
  RNA_target_mouse = DNAStringSet(
    getBM(attributes = c("gene_exon_intron","ensembl_gene_id"),
          filters = "ensembl_gene_id",
          values =
            ortho_ENS$mmusculus_homolog_ensembl_gene,
          mart = martMM)$gene_exon_intron)

  
  # ----------------------------------- milestone 12 --------------------------
  print("milestone 12: ")
  
  # Obtain all human polymorphisms for the RNA target
  if (input$polymorphism_input == TRUE) {
    PMs = getBM(
      attributes = c("minor_allele_freq", "chromosome_start"),
      filters = "ensembl_gene_id",
      values = ensembl_ID,
      mart = martHS
    ) %>%
      as_tibble() %>%
      arrange(chromosome_start, desc(minor_allele_freq)) %>%
      filter(!is.na(minor_allele_freq),
             !duplicated(chromosome_start)) %>%
      rename(chr_start = chromosome_start, PM_freq = minor_allele_freq)
  }

  ##If Ensembl is offline and still want to test -> load in manual test data.
  #PMs <- read.csv("~/PMs.csv")

  
  # ----------------------------------- milestone 13 --------------------------
  print("milestone 13: Get human polymorfisms for RNA target")
  
  # Count Nucleobase sequence occurrences
  
  # Get the sequences
  tr = target_annotation$name
  
  # Count the sequences by making it a table
  replica = table(tr)
  
  # Save it as NoRepeats
  target_annotation$NoRepeats = as.vector(replica[tr])
  
  # ----------------------------------- milestone 14 --------------------------
  print("milestone 14: Count amount of times ASO sequence is repeated in target gene")
  
  # High-Frequency Polymorphisms analysis
  
  # Correcting end and start cord based on direction
  if (chr_coord['strand'] == 1) {
    target_annotation$chr_start = chr_coord['start'] + target_annotation$start - 1
    target_annotation$chr_end = chr_coord['start'] + target_annotation$end - 1
  } else {
    target_annotation$chr_start = chr_coord['end'] - target_annotation$end + 1
    target_annotation$chr_end = chr_coord['end'] - target_annotation$start + 1
  }
  
  # ----------------------------------- milestone 15 -------------------------
  print("milestone 15: Corrected start en end cord based on direction")
  
  # Keep unique names only and extract
  # Information base on chr_start from target.
  if (input$polymorphism_input == TRUE) {
    PM_freq = PMs %>%
      mutate(name = map(chr_start, function(X) {
        filter(target_annotation, chr_start <= X, chr_end >= X) %>%
          select(name, chr_start_anno = chr_start)
      })) %>%
      unnest_legacy() %>%
      rename(chr_start_PM = chr_start) %>%
      group_by(name, chr_start_anno) %>%
      summarise(
        PM_tot_freq = 1 - prod(1 - PM_freq),
        PM_max_freq = max(PM_freq),
        PM_count = n()
      )
    
    # ----------------------------------- milestone 16 -----------------------
    print("milestone 16: Retrieved point mutation frequency count and chance")
    
    # Combine data frames
    target_annotation = left_join(
      target_annotation,
      PM_freq,
      by = c("name" = "name", "chr_start" = "chr_start_anno")
    )
  }
  
  
  # ----------------------------------- milestone 17 -------------------------
  print("milestone 17: Joined PM freq table with target annotation")
  
  # Match RNA Target Regions to the Mouse Ortholog
  if (input$Conserved_input == TRUE) {
    # Get length
    lm = width(RNA_target_mouse)
    
    # Make table of mouse information
    MM_tab = lapply(oligo_lengths, function(i) {
      tibble(st = 1:(lm - i + 1), w = i)
    }) %>%
      bind_rows()
    
    # ----------------------------------- milestone 18.1 -----------------------
    print("milestone 18.1: Match RNA target to mouse ortholog")
    
    # Makes DNAStringSet object with mouse info.
    RNAsitesMM = DNAStringSet(RNA_target_mouse[[1]],
                              start = MM_tab$st,
                              width = MM_tab$w)
    
    # Adds if conserved in mouse.
    target_annotation$conserved_in_mmusculus = target_annotation$name %in% RNAsitesMM
  }

  
  # ----------------------------------- milestone 18.2 -----------------------
  print("milestone 18.2: If selected: Matched mouse ortholog to target gene")
  
  if (input$linux_input == TRUE) {
    # Deze nog toevoegen wanneer viennaRNA werkt
    target_annotation$sec_energy = RNAselffold_R(target_annotation$oligo_seq)
    target_annotation$duplex_energy = RNAduplex_R(target_annotation$oligo_seq)
  }
  # ----------------------------------- milestone 19 -----------------------
  print("milestone 19: Calculated secondary and duplex energy of ASO seq")

  
  ###################count the number of pre-mRNA transcripts with a perfect match 

  #this part will take some time to run...
  uni_tar = dplyr::select(target_annotation, name, length)%>%
    unique() %>%
    split(.,.$length) %>%
    head(10)
  
  uni_tar = lapply(uni_tar, function(X){
    dict0 = PDict(X$name, max.mismatch = 0)
    dict1 = PDict(X$name, max.mismatch = 1)
    
    #perfect match count
    pm = vwhichPDict(
      pdict = dict0, subject = HS,
      max.mismatch = 0, min.mismatch=0)
    X$gene_hits_pm = tabulate(unlist(pm),nbins=nrow(X))
    
    #single mismatch count, without indels
    mm1 = vwhichPDict(
      pdict = dict1, subject = HS,
      max.mismatch = 1, min.mismatch=1)
    X$gene_hits_1mm = tabulate(unlist(mm1),nbins=nrow
                               (X))
    X
  }) %>%
    bind_rows()

  # ----------------------------------- milestone 20 -----------------------
  print("milestone 20: Matched ASO sequences to potential off-targets (perfect match, one mismatch")
  
  target_annotation = left_join(target_annotation, uni_tar, by = c('name', 'length'))
  
  prefilter <- nrow(target_annotation)
  target_annotation <- target_annotation %>%
    filter(
      gene_hits_pm < input$numeric_input_a,
      gene_hits_1mm < input$numeric_input_b
    )

  postfilter <- nrow(target_annotation)
  removed <- prefilter - postfilter

  print(paste0("Filtering Oligo sequences with perfect match > ", input$numeric_input_a, " and 1 mismatch > ", input$numeric_input_b))
  print(paste0("Rows before filtering: ", prefilter))
  print(paste0("Rows after filtering: ", postfilter))
  print(paste0("Filtering removed ", removed, " possible ASOs."))
  
  # ----------------------------------- milestone 21 -----------------------
  print("milestone 21: Filtered ASOs with too many off targets")
  
  # Count the number of pre-mRNA transcripts with a perfect match
  # This part will take some time to run...
  
  summary_server <- target_annotation %>%
    head(2) %>% # For quick off-target testing, use head here
    mutate(results = map2(name, length, ~ {
      res <- all_offt(.x, 2)
      res$name <- .x
      res$length <- .y
      res
    })) %>%
    pull(results) %>%
    bind_rows() %>%
    mutate(distance = mismatches + deletions + insertions,
            gene_name = str_extract(line, "(?<=\\|)[^;]+")
    )
  
  # ----------------------------------- milestone 22 -----------------------
  print("milestone 22: GGGenome searched for all ASO off-targets")
  
  tmp <- tempfile(fileext = ".bgz")
  
  # Download GnomAD lof metrics by gene
  download.file(
    url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
    destfile = tmp,
    mode = "wb"
  )
  gnomad_df <- read_tsv(tmp, show_col_types = FALSE,
                        col_select = c(gene, transcript, oe_lof))
  unlink(tmp)
  
  # ----------------------------------- milestone 23 -----------------------
  print("milestone 23: Downloaded GnomAD lof metrics by gene to dataframe")
  
  off_targets_total = left_join(summary_server, gnomad_df, by = c('gene_name' = 'gene'))

  View(off_targets_total)
  
  # ----------------------------------- milestone 22 -------------------------
  print("milestone 22")
  
  dist_counts <- off_targets_total %>%
    group_by(name, distance) %>%
    summarise(n_distance = n(), .groups = "drop") %>%
    mutate(distance = paste0("n_distance_", distance)) %>%  # bijv. n_distance_0, n_distance_1
    pivot_wider(
      names_from  = distance,
      values_from = n_distance,
      values_fill = 0
    )

  # 2. Per name de laagste oe_lof score bepalen
  oe_lof_min <- off_targets_total %>%
    group_by(name) %>%
    summarise(min_oe_lof = min(oe_lof, na.rm = TRUE), .groups = "drop")

  # 3. Deze twee samenvattingen samenvoegen tot één tabel per name
  off_summary <- dist_counts %>%
    left_join(oe_lof_min, by = "name")

  # 4. Met left_join koppelen aan target_annotation op oligo_seq = name
  target_annotation <- target_annotation %>%
    left_join(off_summary, by = c("name" = "name"))

  View(target_annotation)
  # ----------------------------------- milestone 22 -------------------------
  print("milestone 22")
  
  
  
  target_regions <- target_annotation
  
  # Run the filtering
  if (input$polymorphism_input == TRUE) {
    target_regions <- target_regions %>%
      mutate(across(
        c(PM_max_freq, PM_tot_freq, PM_count),
        ~ replace_na(., 0.0)
      ))
  }
  
  print("temp milestone")
  target_region_select <- target_regions
  
  if (input$perfect_input == TRUE) {
    target_region_select <- filter_function(
      target_region_select,
      input$numeric_input_a,
      "gene_hits_pm",
      input$dropdown_input_a
    )
  }
  if (input$mismatch_input == TRUE) {
    target_region_select <- filter_function(
      target_region_select,
      input$numeric_input_b,
      "gene_hits_1mm",
      input$dropdown_input_b
    )
  }
  if (input$linux_input == TRUE) {
    if (input$Accessibility_input == TRUE) {
      target_region_select <- filter_function(
        target_region_select,
        input$numeric_input_c,
        "accessibility",
        input$dropdown_input_c
      )
    }
  }
  if (input$polymorphism_input == TRUE) {
    if (input$Poly_input == TRUE) {
      target_region_select <- filter_function(
        target_region_select,
        input$numeric_input_d,
        "PM_tot_freq",
        input$dropdown_input_d
      )
    }
  }
  if (input$tox_input == TRUE) {
    target_region_select <- filter_function(
      target_region_select,
      input$numeric_input_e,
      "tox_score",
      input$dropdown_input_e
    )
  }
  if (input$Conserved_input == TRUE) {
    target_region_select <- filter(target_region_select, conserved_in_mmusculus == TRUE)
  }
  
  # Check if the filtered result is empty
  if (nrow(target_region_select) == 0) {
    # If empty, return the original target_regions
    target_region_select <- target_regions
    print("No results after filtering, exporting unfiltered list")
    showNotification(
      "No results after filtering, exporting unfiltered list",
      type = "default",
      duration = NULL,
      # Notification stays until clicked away
      closeButton = TRUE
    ) # Include a close button)
  }
  
  # ----------------------------------- milestone 19 -------------------------
  print("milestone19")
  
  # Get all unique starting positions in order
  start_pos = sort(unique(target_region_select$start))
  
  num_data_points <- length(start_pos)
  
  # ----------------------------------- milestone 19.1 -----------------------
  print("milestone19.1")
  
  # Set amount of clusters
  if (num_data_points > 1) {
    K <- min(10, num_data_points - 1)
  } else {
    # Handle the case when there are not enough data points
    K <- 1
  }
  
  # ----------------------------------- milestone 19.2 -----------------------
  print("milestone19.2")
  
  # Cluster the selected regions
  cluster_tab = tibble(
    start = start_pos,
    cluster = clara(
      x = start_pos,
      k = K,
      metric = 'euclidean',
      pamLike = T,
      samples = 100
    )$clustering
  )
  
  # ----------------------------------- milestone 20 -------------------------
  print("milestone20")
  
  nucleobase_select = left_join(target_region_select, cluster_tab, by =
                                  'start') %>%
    group_by(cluster) %>%
    sample_n(1) %>%
    ungroup()
  
  # ----------------------------------- milestone 21 -------------------------
  print("milestone21")
  
  # Render the tables.
  output$results1 <- renderDataTable({
    datatable(target_region_select,
              rownames = FALSE,
              selection = "single") %>%
      formatStyle(names(target_region_select), color = "black")
  })
  
  output$results2 <- renderDataTable({
    datatable(nucleobase_select,
              rownames = FALSE,
              selection = "single") %>%
      formatStyle(names(nucleobase_select), color = "black")
  })
  
  # ----------------------------------- Fayes deel ---------------------------
  
  current_seq <- reactiveVal(NULL)
  current_mismatch <- reactiveVal(2)
  current_offtargets <- reactiveVal(NULL)
  cached_results <- reactiveVal(list())
  
  compute_offtarget_accessibility <- function(df) {
    if (input$linux_input != TRUE) return(df)
    
    for (i in seq_len(nrow(df))) {
      offtarget_seq <- gsub("-", "", df$subject_seq[i])
      l_ot <- nchar(offtarget_seq)
      
      snip_result <- acc_snippet(
        begin_target  = df$start_target[i],
        end_target    = df$end_target[i],
        begin_snippet = df$snippet_start[i],
        end_snippet   = df$snippet_end[i],
        full_snippet  = df$snippet[i]
      )
      
      l_snipseq <- nchar(snip_result$snippet_seq)
      
      df$offtarget_accessibility[i] <-
        RNAplfold_R(
          snip_result$snippet_seq,
          u.in = l_ot
        ) %>%
        as_tibble() %>%
        mutate(end = 1:l_snipseq) %>%
        gather(length, accessibility, -end) %>%
        mutate(length = as.integer(length)) %>%
        filter(
          length == l_ot,
          (end - l_ot + 1) <= snip_result$target_end_internal,
          end >= snip_result$target_start_internal
        ) %>%
        summarise(mean_accessibility = mean(accessibility, na.rm = TRUE)) %>%
        pull(mean_accessibility)
    }
    
    return(df)
  }
  
      observeEvent(input$apply_mismatch, {
        req(current_seq())
        
        seq <- toupper(current_seq())
        mm  <- as.numeric(input$user_mismatch)
        key <- paste0("mm", mm)
        
        current_mismatch(mm)
        cache <- cached_results()
        
        if (!is.null(cache[[seq]]) && !is.null(cache[[seq]][[key]])) {
          subset_df <- cache[[seq]][[key]]
          
        } else {
          if (mm %in% c(0, 1, 2)) {
            subset_df <- summary_server %>%
              filter(toupper(name) == seq, mismatches <= mm)
            
          } else if (mm == 3) {
            showNotification("Results loading", type = "default",
                             duration = NULL, closeButton = TRUE)
            new_res <- all_offt(seq, 3)
            print("mm3 gedaan")
            new_res$name   <- seq
            new_res$length <- nchar(seq)
            subset_df <- new_res
          }
          
          subset_df <- compute_offtarget_accessibility(subset_df)
          
          if (is.null(cache[[seq]])) cache[[seq]] <- list()
          cache[[seq]][[key]] <- subset_df
          cached_results(cache)
        }
        
        if (input$linux_input == TRUE) {
          subset_df <- subset_df %>% 
            relocate(offtarget_accessibility, .after = insertions)
        }
        
        current_offtargets(subset_df)
        
        output$numb_offtargets <- renderText(
          paste0("# off targets: ", nrow(subset_df))
        )
      })
      
      output$offtarget_results <- DT::renderDataTable({
        req(current_offtargets())
        datatable(current_offtargets(), rownames = FALSE)
      })
      
      output$download_offtarget <- downloadHandler(
        filename = function() {
          paste(
            'offtargets_', current_seq(), "_mismatches_", current_mismatch(),
            "_", Sys.Date(), '.csv'
          )
        },
        content = function(con) {
          write.csv2(current_offtargets(), con, row.names = FALSE)
        }
      )
  
  # ----------------------------------- Harrys deel --------------------------
  
  # String reverse function
  reverse_string <- function(x) {
    paste0(rev(strsplit(x, "")[[1]]), collapse = "")
  }
  
  # A stored value for use in the second table and download.
  selected_target <- reactiveVal(NULL)
  oligo_sequence <- reactiveVal(NULL)
  rnaseh_stored <- reactiveVal(NULL)
  
  # Apply end modifications.
  observeEvent(input$add_mods, {
    row_data <- selected_target()
    if (is.null(row_data))
      return()
    
    rnaseh_data <- rnaseh_results(
      selected_row_name = row_data$name,
      oligo_seq = row_data$oligo_seq,
      mod_5prime = input$mod_5prime,
      mod_3prime = input$mod_3prime
    )
    
    rnaseh_stored(rnaseh_data)
    
    output$rnaseh_results <- renderDataTable({
      datatable(rnaseh_data, selection = list(mode = 'single', selected = 1))
    })
    
    output$cleavage_visual <- renderUI(div())
  })
  
  # The second observer object gives a visual of the cleavage site on the target sequence.
  observeEvent(input$rnaseh_results_rows_selected, {
    row_number <- input$rnaseh_results_rows_selected
    if (length(row_number) == 0)
      return()
    
    row_data <- selected_target()
    if (is.null(row_data))
      return()
    
    rnaseh_data <- rnaseh_stored()
    if (is.null(rnaseh_data))
      return()
    
    selected_row <- rnaseh_data[row_number, ]
    
    # Oligo sequence
    oligo_seq <- row_data$oligo_seq
    
    mod5 <- input$mod_5prime
    mod3 <- input$mod_3prime
    
    oligo_len <- nchar(oligo_seq)
    
    mod5_region <- substr(oligo_seq, 1, mod5)
    mod3_region <- substr(oligo_seq, oligo_len - mod3 + 1, oligo_len)
    
    mid_region <- substr(oligo_seq, mod5 + 1, oligo_len - mod3)
    
    # RNA sequence
    position_string <- str_split_1(selected_row$position, " - ")
    start_pos <- as.numeric(position_string[1])
    
    target_seq_fw <- row_data$name
    target_seq_rv <- reverse_string(target_seq_fw)
    
    rna_len <- nchar(target_seq_fw)
    
    cleavage_start_fw <- start_pos
    cleavage_pos_fw <- cleavage_start_fw + 6
    cleavage_pos_rv <- rna_len - cleavage_pos_fw
    cleavage_site_up <- cleavage_pos_rv - 1
    cleavage_site_down <- cleavage_pos_rv + 7
    
    upstream <- substr(target_seq_rv, 1, cleavage_site_up - 1)
    site_start <- substr(target_seq_rv, cleavage_site_up, cleavage_pos_rv - 1)
    cut_site <- substr(target_seq_rv, cleavage_pos_rv, cleavage_pos_rv)
    site_down <- substr(target_seq_rv, cleavage_pos_rv + 1, cleavage_site_down)
    downstream <- substr(target_seq_rv, cleavage_site_down + 1, rna_len)
    
    oligo_visual_fw <- paste0(
      "<b style='color:darkorange;'>5'</b> ",
      "<span style='font-weight:bold; color:#90D5FF;'>",
      mod5_region,
      "</span>",
      mid_region,
      "<span style='font-weight:bold; color:#90D5FF;'>",
      mod3_region,
      "</span>",
      " <b style='color:darkorange;'>3'</b>"
    )
    
    rna_visual_rv <- paste0(
      "<b style='color:darkorange;'>3'</b> ",
      upstream,
      "<span style='background-color: lightblue; color: red; font-weight: bold;'>",
      site_start,
      cut_site,
      "|",
      site_down,
      "</span>",
      downstream,
      " <b style='color:darkorange;'>5'</b>"
    )
    
    output$cleavage_visual <- renderUI({
      HTML(
        paste0(
          "<h5>Oligo sequence (ASO): </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 18px;'>",
          oligo_visual_fw,
          "</div>",
          "<h5>RNA sequence: </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 18px;'>",
          rna_visual_rv,
          "</div>"
        )
      )
    })
  })
  
  # Download handler.
  output$download_rnaseh <- downloadHandler(
    filename = function() {
      row_data <- selected_target()
      if (is.null(row_data)) {
        "RNaseH_results.xlsx"
      } else {
        paste0("RNaseH_results_", row_data$name, ".xlsx")
      }
    },
    content = function(file) {
      data <- rnaseh_stored()
      
      if (is.null(data))
        data <- data.frame(Message = "No avalable data")
      
      data[] <- lapply(data, as.character)
      
      wb <- createWorkbook()
      addWorksheet(wb, "RNaseH_results")
      writeData(wb, "RNaseH_results", data)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ----------------------------------- Table handlers -----------------------
  
  # Function for table handler call
  handle_table_events <- function(input,
                                  output,
                                  session,
                                  table_id,
                                  table_data,
                                  off_targets_total,
                                  selected_target,
                                  oligo_sequence,
                                  rnaseh_stored,
                                  current_seq,
                                  current_mismatch,
                                  current_offtargets,
                                  cached_results) {
    observeEvent({
      list(input[[paste0(table_id, "_cell_clicked")]], input[[paste0(table_id, "_rows_selected")]])
    }, {
      cell <- input[[paste0(table_id, "_cell_clicked")]]
      row  <- input[[paste0(table_id, "_rows_selected")]]
      
      if (!is.null(cell) && !is.null(cell$row)) {
        session$sendCustomMessage("selectRow", list(table = table_id, row   = cell$row))
      }
      
      if (!is.null(row) && length(row) > 0) {
        other_table <- if (table_id == "results1")
          "results2"
        else
          "results1"
        proxy_other <- dataTableProxy(other_table)
        selectRows(proxy_other, NULL)
        
        row_data <- table_data[row, ]
        
        # Off-target functionality
              req(row_data$name)
              seq <- toupper(row_data$name)
              
              if (!grepl("^[ACGT]+$", seq)) return()
              current_seq(seq)
              mm <- 2
              current_mismatch(mm)
              updateSelectInput(session, "user_mismatch", selected = 2)
              
              key <- paste0("mm", mm)
              cache <- cached_results()
              
              if (!is.null(cache[[seq]]) && !is.null(cache[[seq]][[key]])) {
                default_subset <- cache[[seq]][[key]]
                
              } else {
                default_subset <- summary_server %>%
                  filter(toupper(name) == seq, mismatches <= mm)
              }
              
              if (input$linux_input == TRUE) {
                default_subset <- compute_offtarget_accessibility(default_subset)
                default_subset <- default_subset %>% 
                  relocate(offtarget_accessibility, .after = insertions)
              }
              
              # Cache
              cache[[seq]][[key]] <- default_subset
              cached_results(cache)
              
              # results
              current_offtargets(default_subset)
              
              output$offtarget_title <- renderText(
                paste0("Off target results for: ", seq)
              )
              output$aso_seq <- renderText(
                paste0("ASO sequence: ", as.character(reverseComplement(DNAString(seq))))
              )
              output$numb_offtargets <- renderText(
                paste0("# off targets: ", nrow(default_subset))
              )
       
        # RNaseH functionality
        selected_target(row_data)
        oligo_sequence(row_data$oligo_seq)
        
        rnaseh_data <- rnaseh_results(
          selected_row_name = row_data$name,
          mod_5prime = 0,
          mod_3prime = 0
        )
        
        rnaseh_stored(rnaseh_data)
        
        output$rnaseh_title <- renderText(paste0("RNase H results for: ", row_data$name))
        
        output$rnaseh_info <- renderText({
          HTML(
            paste0(
              "length of sequence: ",
              row_data$length,
              "<br>",
              "Oligo sequence (ASO): ",
              oligo_sequence()
            )
          )
        })
        
        output$rnaseh_results <- renderDataTable({
          datatable(rnaseh_data,
                    selection = list(mode = "single", selected = 1))
        })
        
        updateTabsetPanel(session, "tabs_main", selected = "RNase H cleavage results")
      }
      
    })
  }
  
  # Table one handler call
  handle_table_events(
    input = input,
    output = output,
    session = session,
    table_id = "results1",
    table_data = target_region_select,
    off_targets_total = off_targets_total,
    selected_target = selected_target,
    oligo_sequence = oligo_sequence,
    rnaseh_stored = rnaseh_stored,
    current_seq = current_seq,
    current_mismatch = current_mismatch,
    current_offtargets = current_offtargets,
    cached_results = cached_results
  )
  
  # Table two handler call
  handle_table_events(
    input = input,
    output = output,
    session = session,
    table_id = "results2",
    table_data = nucleobase_select,
    off_targets_total = off_targets_total,
    selected_target = selected_target,
    oligo_sequence = oligo_sequence,
    rnaseh_stored = rnaseh_stored,
    current_seq = current_seq,
    current_mismatch = current_mismatch,
    current_offtargets = current_offtargets,
    cached_results = cached_results
  )
  
  # ----------------------------------- End of script ------------------------
  
   # Collect the data, change the name for each gene tested
      output$Download_input <- downloadHandler(

        filename = function() {
          current_date <- format(Sys.time(), "%Y-%m-%d %H-%M-%S")
          paste0("Result_output_", current_date, ".zip")
        },

        content = function(file) {
          file1 <- tempfile(fileext = ".csv")
          file2 <- tempfile(fileext = ".csv")

          write.csv(target_region_select(), file1, row.names = FALSE)
          write.csv(nucleobase_select(), file2, row.names = FALSE)

          zip(file, c(file1, file2), flags = "-j")
        }
      )

  t2 <- Sys.time()
  time <- t2 - t1
  print(time)
  showNotification(
    "Script finished",
    type = "default",
    duration = NULL,
    # Notification stays until clicked away
    closeButton = TRUE
  )
  print("done")
  }
  })
}


