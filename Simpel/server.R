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
library(DT)
library(shinyBS)
#setwd("/home/gebruiker")
# setwd("..")

setwd("C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2")
print("work directory set")

source("tools/offtarget.R")

function(input, output, session) {
  showNotification("Finished loading", type = "default", duration = NULL, # Notification stays until clicked away
                   closeButton = TRUE) # Include a close button)
  observeEvent(input$run_button, {
    showNotification("Script started", type = "default", duration = NULL, # Notification stays until clicked away
                     closeButton = TRUE) # Include a close button)
    #######################################################################functions
    #Make the tox score function
    calculate_acute_neurotox <- function(xx) {
      ## make sure input is in character format
      xx <- as.character(xx)
      
      ## count the number of each nucleotide
      lf <- function(x) {
        x <- tolower(x)
        x <- strsplit(x, "")[[1]]
        x <- table(factor(x,levels=c("a","c","t","g")))
        return(x)
      }
      cnt_nt <- as.data.frame(t(sapply(xx, lf)))
      ## count number of nucleotides from the 3'-end untill it finds the first g
      gfree3 <- function(x) {
        x <- tolower(x)
        x <- strsplit(x, "")[[1]]
        tfg <- x=="g"
        if (sum(tfg)==0) {
          l3 <- NA
        } else {
          posg <- c(1:length(x))[tfg]
          l3 <- length(x)-max(posg)
        }
        return(l3)
      }
      cnt_gfree3 <- sapply(xx, gfree3)
      cnt_gfree3[cnt_gfree3>20] <- 20      #Set max to 20
      cnt_gfree3[is.na(cnt_gfree3)] <- 20  #Set no g in ASO to 20
      ## Calculate final score based on trained parameters and return result
      calc_out <- round(136.0430 - 3.1263*cnt_nt$a - 5.1100*cnt_nt$c -
                          4.7217*cnt_nt$t - 10.1264*cnt_nt$g + 1.3577*cnt_gfree3,1) 
      
      return(as.numeric(calc_out))
    }
    
    if (input$linux_input == TRUE) {
    #2.5.1 Predict accessibility for an RNA target
    RNAplfold_R = function(seq.char, L.in = 40,W.in = 80, u.in = 16){
      cmmnd2 = paste("RNAplfold -L", L.in, "-W", W.in, "-u",u.in)
      seq.char = as.character(seq.char)
      cat(seq.char, file = paste("|", cmmnd2, sep = ""))
      acc.tx = read.delim("plfold_lunp", as.is = T, skip = 2,header = F, row.names = 1)
      acc.tx = acc.tx[, colSums(is.na(acc.tx)) !=nrow(acc.tx)]
      colnames(acc.tx) = 1:ncol(acc.tx)
      file.remove("plfold_lunp")
      file.remove("plfold_dp.ps")
      return(acc.tx)
    }
    
    #2.5.3 Predict duplex formation and self folding of oligonucleotides
    RNAduplex_R = function(seqs){
      sys_cmd = system('RNAduplex',input = c(seqs,seqs),intern = TRUE)
      as.numeric(regmatches(sys_cmd, regexpr("-?\\d+\\.\\d+", sys_cmd)))
    }
    
    RNAselffold_R = function (seqs) {
      output = system('RNAfold --noPS', input=c(seqs),intern = T)
      output = unlist(strsplit(output[grepl('[0-9]', output)],'[(]'))
      as.double(gsub(' |[)]', '', output[grepl('[0-9]',output)]))
    }
    }
    
    filter_function <- function(df, valueX, columname, operator_string) {
      switch(operator_string,
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
    ##############################################Store all human pre-mRNA sequences
    # path = getwd()
    # getwd()
    
    path <- "C:/Users/fayef/Documents/BI/BI3/periode_1/XEXT/ASOstool-v2"
    # Load the TxDb object
    txdb_hsa <- loadDb("txdb_hsa_biomart.db")
    
    print("milestone1")
    
    # Extract the genes
    gdb_hsa <- genes(txdb_hsa)
    
    # Define the chromosomes to keep
    chr_to_keep <- c(as.character(1:22),'X','Y','MT')
    
    # filtert Hsapiens zodat het alleen de aangegeven chromosomen pakt
    Hsapiens@user_seqnames <- setNames(chr_to_keep,chr_to_keep)
    Hsapiens@seqinfo <- Hsapiens@seqinfo[chr_to_keep]
    
    # Subset genes to keep only those on specified chromosomes
    gdb_hsa <- gdb_hsa[seqnames(gdb_hsa) %in% chr_to_keep]
    
    print("milestone2")
    
    # Get the sequences *
    HS <- getSeq(Hsapiens, gdb_hsa)
    
    ############################################target collect the pre-mRNA sequence
    
    #define wanted Ensembl ID
    print("milestone3")
    ensembl_ID = input$ensemble_id_input
    
    #retrieve a specific RNA target using the Ensembl ID
    RNA_target = HS[names(HS)==ensembl_ID]
    print("milestone4")
    
    #filters on ensembl ID
    target_ranges = gdb_hsa[names(gdb_hsa)==ensembl_ID]
    
    #extracts the chromosome name for target range,extracts the start and end 
    #position of the genomic region and note from which strand it is. 
    #positive strand 1 ('+'), negative strand -1 ('-'), or unspecified 0 ('*').
    chr_coord = c(
      chr=as.numeric(as.character(seqnames(target_ranges))),
      start=start(target_ranges),
      end=end(target_ranges),
      strand=ifelse(strand(target_ranges)=="+",1,-1))
    
    
    print("milestone5")

    ###############################Obtain the mouse ortholog of the human RNA target
    
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
    
    print("milestone6")
    
    RNA_target_mouse = DNAStringSet(
      getBM(attributes = c("gene_exon_intron","ensembl_gene_id"),
            filters = "ensembl_gene_id",
            values =
              ortho_ENS$mmusculus_homolog_ensembl_gene,
            mart = martMM)$gene_exon_intron)
    }
    
    print("milestone7")
    
    ###############################Obtain all human polymorphisms for the RNA target
    if (input$polymorphism_input == TRUE) {
      PMs = getBM(attributes = c("minor_allele_freq",
                               "chromosome_start"),
                filters = "ensembl_gene_id",
               values = ensembl_ID,
                mart = martHS) %>%
      as_tibble() %>%
      arrange(chromosome_start, desc(minor_allele_freq)) %>%
      filter(!is.na(minor_allele_freq), !duplicated(chromosome_start)) %>%
      rename(chr_start=chromosome_start, PM_freq=minor_allele_freq)
    }
    
    ##If Ensembl is offline and still want to test -> load in manual test data.
    #PMs <- read.csv("~/PMs.csv")
    
    print("milestone8")
    
    #set the width of rna_target as value L 
    l = width(RNA_target)
    
    #note length of subsequence
    oligo_lengths = input$oligo_length_range[1]:input$oligo_length_range[2]
    
    #construct the dataframe
    target_annotation = lapply(oligo_lengths, function(i){
      tibble(start=1:(l-i+1),length=i)}) %>%
      bind_rows() %>%
      mutate(end = start+length-1)

    #Adds "name" sequence to the dataframe
    target_regions = DNAStringSet(
      RNA_target[[1]],
      start=target_annotation$start,
      width=target_annotation$length)
    names(target_regions) = as.character(target_regions)
    
    target_annotation$name = names(target_regions)
    
    print("milestone9")
    
    ###########################################Count Nucleobase sequence occurrences
    
    #get the sequences
    tr = target_annotation$name
    
    #count the sequences by making it a table
    replica = table(tr)
    
    #save it as NoRepeats
    target_annotation$NoRepeats = as.vector(replica[tr])
    
    
    print("milestone10")
    
    ###################count the number of pre-mRNA transcripts with a perfect match 
    
    #this part will take some time to run...
    
    target_annotation <- target_annotation %>%
      filter(!grepl("^C", name)) %>%
      mutate(reverse_comp = reverse_complement(name))

    summary_server <- target_annotation %>%
      head(2) %>%
      mutate(results = map2(name, length, ~ {
        res <- all_offt(.x, 2)         
        res$name <- .x          
        res$length <- .y        
        res                    
      })) %>%
      pull(results) %>%
      bind_rows()
    
    
    print(summary_server)

    uni_tar <- summary_server %>%
      group_by(name, length) %>%
      summarise(
        gene_hits_pm  = sum(`Total Mismatches` == 0, na.rm = TRUE),
        gene_hits_1mm = sum(`Total Mismatches` == 1, na.rm = TRUE),
        .groups = "drop"
      )

    
    print("milestone11")
    
    if (input$linux_input == TRUE) {
    #3.4 Estimate Transcript Accessibility for the RNA Target at Single-Nucleotide Resolution
    accessibility = RNAplfold_R(RNA_target,u.in = max(oligo_lengths)) %>%
      as_tibble() %>%
      mutate(end=1:l) %>%
      gather(length,accessibility, -end) %>%
      mutate(length=as.double(length))
    
    target_annotation = left_join(target_annotation,accessibility,by=c('length','end'))
    }
    print("milestone11.5")
    ###########################################High-Frequency Polymorphisms analysis
    
    #correcting end and start cord based on direction   
    if (chr_coord['strand'] == 1) {
      target_annotation$chr_start = chr_coord['start'] + target_annotation$start - 1
      target_annotation$chr_end = chr_coord['start'] + target_annotation$end - 1
    } else {
      target_annotation$chr_start = chr_coord['end'] - target_annotation$end + 1
      target_annotation$chr_end = chr_coord['end'] - target_annotation$start + 1
    }
    
    print("milestone12")
    
    #keep unique names only and extract 
    #information base on chr_start from target.
    if (input$polymorphism_input == TRUE) {
    PMmax = PMs %>%
      mutate(name = map(chr_start, function(X){
        filter(target_annotation,
               chr_start <= X,
               chr_end >= X)$name
      })) %>%
      unnest_legacy() %>%
      arrange(desc(PM_freq)) %>%
      filter(!duplicated(name))
    
    print("milestone13")
    
    #combine data frames
    target_annotation = left_join(
      target_annotation, PMmax,
      by = c("name", "chr_start"))
    }
    
    print("milestone14")
    
    ##################################Match RNA Target Regions to the Mouse Ortholog
    if (input$Conserved_input == TRUE) {
    #get length
    lm = width(RNA_target_mouse)
    
    #make table of mouse information
    MM_tab = lapply(oligo_lengths, function(i){
      tibble(st=1:(lm-i+1),w=i)}) %>%
      bind_rows()
    
    print("milestone15")
    
    #makes dnastringset object with mouse info.
    RNAsitesMM = DNAStringSet(
      RNA_target_mouse[[1]],
      start=MM_tab$st,
      width=MM_tab$w)
    
    #adds if conserved in mouse.
    uni_tar$conserved_in_mmusculus = uni_tar$name %in% RNAsitesMM
    }
    print("milestone16")
    
    #bereken aantal "cg"
    uni_tar$CGs = ( uni_tar$length -
                      nchar(gsub('CG','',uni_tar$name)) )/2
    
    #maak de unique reverse complements van sequences in target regions
    nucleobase_seq =unique(reverseComplement(target_regions))
    
    #voeg ze toe aan uni_tar gekoppeld aan name
    uni_tar$oligo_seq = as.character(nucleobase_seq[uni_tar$name])
    
    #toxicity score acute neurotox score (desired >70)
    uni_tar$tox_score = calculate_acute_neurotox(uni_tar$oligo_seq)
    
    print("milestone17.1")
    
    if (input$linux_input == TRUE) {
    #### deze nog toevoegen wanneer viennaRNA werkt
    uni_tar$sec_energy = RNAselffold_R(uni_tar$oligo_seq)
    uni_tar$duplex_energy = RNAduplex_R(uni_tar$oligo_seq)
  
    print("milestone17.2")
    }
    #join tables
    target_regions = left_join(
      target_annotation,
      uni_tar,
      by=c('name','length'))
    
    print("milestone18")
    
    #run the filtering
     if (input$polymorphism_input == TRUE) {
      target_regions <- target_regions %>%
      mutate(PM_freq  = replace_na(PM_freq , 0.0))
    }
    
    print("temp milestone")
    target_region_select <- target_regions
    
    if (input$perfect_input == TRUE) {
      target_region_select <- filter_function(target_region_select, input$numeric_input_a, "gene_hits_pm", input$dropdown_input_a)
    }
    if (input$mismatch_input == TRUE) {
      target_region_select <- filter_function(target_region_select, input$numeric_input_b, "gene_hits_1mm", input$dropdown_input_b)
    }
    if (input$linux_input == TRUE) {
      if (input$Accessibility_input == TRUE) {
      target_region_select <- filter_function(target_region_select, input$numeric_input_c, "accessibility", input$dropdown_input_c)
      }
    }
    if (input$polymorphism_input == TRUE) {
     if (input$Poly_input == TRUE) {
        target_region_select <- filter_function(target_region_select, input$numeric_input_d, "PM_freq", input$dropdown_input_d)
      }
    }
    if (input$tox_input == TRUE) {
      target_region_select <- filter_function(target_region_select, input$numeric_input_e, "tox_score", input$dropdown_input_e)
    }
    if (input$Conserved_input == TRUE) {
      target_region_select <- filter(target_region_select, conserved_in_mmusculus == TRUE)
    }
    
    
    # Check if the filtered result is empty
    if (nrow(target_region_select) == 0) {
      # If empty, return the original target_regions
      target_region_select <- target_regions
      print("No results after filtering, exporting unfiltered list")
      showNotification("No results after filtering, exporting unfiltered list", type = "default", duration = NULL, # Notification stays until clicked away
                       closeButton = TRUE) # Include a close button)
    }
    
    ############################################################################

    
    print("milestone19")
    
    #get all unique starting positions in order
    start_pos = sort(unique(target_region_select$start))
    
    num_data_points <- length(start_pos)
    
    print("milestone19.1")
    
    #Set amount of clusters
    if (num_data_points > 1) {
      K <- min(10, num_data_points - 1)
    } else {
      # Handle the case when there are not enough data points
      K <- 1
    }
    
    print("milestone19.2")
    
    
    #cluster the selected regions
    cluster_tab = tibble(
      start = start_pos,
      cluster = clara(x = start_pos, k = K,
                      metric = 'euclidean',
                      pamLike = T, samples=100)$clustering)
    
    print("milestone20")
    
    nucleobase_select = left_join(target_region_select,
                                  cluster_tab, by='start') %>%
      group_by(cluster) %>%
      sample_n(1) %>%
      ungroup()
    
    print("milestone21")
    
    output$results1 <- renderDataTable({
      datatable(target_region_select, rownames = FALSE) %>%
        formatStyle(names(target_region_select), color = "black")
    })

    current_seq <- reactiveVal(NULL)
    current_mismatch <- reactiveVal(2)
    current_offtargets <- reactiveVal(NULL)
    cached_results <- reactiveVal(list())
    
    observeEvent(input$results1_cell_clicked, {
      
      seq <- toupper(input$results1_cell_clicked$value)
      req(seq)
      req(grepl("^[ACGT]+$", seq, ignore.case = TRUE))
      
      current_seq(seq)
      current_mismatch(2)
      updateSelectInput(session, "user_mismatch", selected = 2)
      
      default_subset <- summary_server %>%
        filter(toupper(name) == seq, `Total Mismatches` <= 2)
      
      current_offtargets(default_subset)
      
      output$offtarget_title <- renderText(paste0("Off target results for: ", seq))
      output$aso_seq <- renderText(paste0("ASO sequence: ", as.character(reverseComplement(DNAString(seq)))))
      output$numb_offtargets <- renderText(paste0("# off targets: ", nrow(default_subset)))
    })
    
    observeEvent(input$apply_mismatch, {
      
      req(current_seq())
      
      mm  <- as.numeric(input$user_mismatch)
      seq <- toupper(current_seq())
      key <- paste0("mm", mm)
      
      current_mismatch(mm)
      cache <- cached_results()
      
      if (!is.null(cache[[seq]]) && !is.null(cache[[seq]][[key]])) {
        subset_df <- cache[[seq]][[key]]
        
      } else {
        
        
        if (mm %in% c(0,1,2)) {
          
          subset_df <- summary_server %>%
            filter(toupper(name) == seq, `Total Mismatches` <= mm)
          
        } else if (mm == 3) {
          showNotification("Results loading", type = "default", duration = NULL, # Notification stays until clicked away
                           closeButton = TRUE)
          
          new_res <- all_offt(seq, 3)
          new_res$name   <- seq
          new_res$length <- nchar(seq)
          
          subset_df <- new_res
        }
        
        if (is.null(cache[[seq]])) cache[[seq]] <- list()
        
        cache[[seq]][[key]] <- subset_df
        
        cached_results(cache)
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
        paste('offtargets_', current_seq(), "_mismatches_", current_mismatch(), "_", Sys.Date(), '.csv')
      },
      content = function(con) {
        data_offtarget <- current_offtargets()
        req(data_offtarget)
        write.csv2(data_offtarget, con, row.names = FALSE)
      }
    )
    
    # -----------------------------------

    
    output$results2 <- renderDataTable({
      datatable(nucleobase_select, rownames = FALSE) %>%
        formatStyle(names(nucleobase_select), color = "black")
    })
    
    # Get the current date
    current_date <- format(Sys.time(), "%Y-%m-%d %H-%M-%S")
    
    
    #collect the data, change the name for each gene tested
    if (input$Download_input == TRUE) {
    write.csv(target_region_select, file = paste0("Results_output_", current_date, ".csv"))
    write.csv(nucleobase_select, file = paste0("Results_output_clustered_", current_date, ".csv"))
    }
    showNotification("Script finished", type = "default", duration = NULL, # Notification stays until clicked away
                     closeButton = TRUE) # Include a close button)
    print("done")
  })
  
}