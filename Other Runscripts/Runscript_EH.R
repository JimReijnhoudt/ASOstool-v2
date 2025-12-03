##############################################
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(GenomicFeatures)
library(AnnotationDbi)
library(BiocManager)
library(cluster)

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
  ## count number of nucleotides from the 3'-end that are not g
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

############################################################################path
path = "/mnt/data/Jeremy/txdb/"

# Load the TxDb object
txdb_hsa <- loadDb("/opt/ASOstool-v2/txdb_hsa_biomart.db")

# Extract the genes *
gdb_hsa <- genes(txdb_hsa)

# Define the chromosomes to keep
chr_to_keep <- c(as.character(1:22),'X','Y','MT')  ####

# filtert Hsapiens zodat het alleen de aangegeven chromosomen pakt
Hsapiens@user_seqnames <- setNames(chr_to_keep,chr_to_keep)
Hsapiens@seqinfo <- Hsapiens@seqinfo[chr_to_keep]

# Subset genes to keep only those on specified chromosomes
gdb_hsa <- gdb_hsa[seqnames(gdb_hsa) %in% chr_to_keep]

# Get the sequences *
HS <- getSeq(Hsapiens, gdb_hsa)

############################################target collect the pre-mRNA sequence
#define wanted Ensembl ID

#ensembl_ID = "SCN2A"
ensembl_ID = 'ENSG00000136531'

#retrieve a specific RNA target using the Ensembl ID
RNA_target = HS[names(HS)==ensembl_ID]

#############preparation polymorphisms analysis organize chromosomal coordinates

#filters on ensembl ID
target_ranges = gdb_hsa[names(gdb_hsa)==ensembl_ID]


#extracts the chromosome name for target range,extracts the start and end 
#position of the genomic region and note from which strand it is. positive strand 1 ('+'), negative strand -1 ('-'), or unspecified 0 ('*').
chr_coord = c(
  chr=as.numeric(as.character(seqnames(target_ranges))),
  start=start(target_ranges),
  end=end(target_ranges),
  strand=ifelse(strand(target_ranges)=="+",1,-1))

###############################Obtain the mouse ortholog of the human RNA target

library(biomaRt)
library(Biostrings)


# Define the marts for mmusculus and hsapiens
martMM = useMart(biomart="ensembl",
                 dataset="mmusculus_gene_ensembl")
martHS = useMart(biomart="ensembl",
                 dataset="hsapiens_gene_ensembl")

# Get the orthologous Ensembl gene for the provided human Ensembl ID
ortho_ENS = getBM(attributes = "mmusculus_homolog_ensembl_gene",
                  filters = "ensembl_gene_id",
                  values = ensembl_ID, mart = martHS,
                  bmHeader = FALSE)

##if "Error in gzfile" occurs, check if the folder exist, else make it.

# Extract the Ensembl IDs of the mouse orthologs and get gene exon/intron information for mouse orthologs
RNA_target_mouse = DNAStringSet(
  getBM(attributes = c("gene_exon_intron","ensembl_gene_id"),
        filters = "ensembl_gene_id",
        values =
          ortho_ENS$mmusculus_homolog_ensembl_gene,
        mart = martMM)$gene_exon_intron)

###############################Obtain all human polymorphisms for the RNA target
library(tidyverse)

#Arranges the data by chromosome_start in ascending order 
#and then by minor_allele_freq in descending order, filters NA and removes dupes
#renames chr_start en pm_freq.

PMs = getBM(attributes = c("minor_allele_freq",
                           "chromosome_start"),
            filters = "ensembl_gene_id",
            values = ensembl_ID,
            mart = martHS) %>%
  as_tibble() %>%
  arrange(chromosome_start, desc(minor_allele_freq)) %>%
  filter(!is.na(minor_allele_freq), !duplicated(chromosome_start)) %>%
  rename(chr_start=chromosome_start, PM_freq=minor_allele_freq)

########
#plan to add viennaRNA on this if linux :D

########
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

plot.postscript.file <- function(file = "Rplots.ps")
{
  # define viewer for UNIX/LINUX or Windows
  viewer <- ifelse(.Platform$OS.type == "unix", "gv", "GSview")
  
  system(paste(viewer, file, sep = " "))
}
#####################################################Construct target_annotation

#set the width of rna_target as value L 
l = width(RNA_target)

#note length of subsequence
oligo_lengths = 14:20    ###


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


###########################################Count Nucleobase sequence occurrences

#get the sequences
tr = target_annotation$name

#count the sequences by making it a table
replica = table(tr)

#save it as NoRepeats
target_annotation$NoRepeats = as.vector(replica[tr])

###################count the number of pre-mRNA transcripts with a perfect match 
start.time1 <- Sys.time()

#this part will take some time to run...
uni_tar = dplyr::select(target_annotation, name, length)%>%
  unique() %>%
  split(.,.$length)

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

end.time1 <- Sys.time()
time.taken <- round(end.time1 - start.time1,2)
time.taken
############################################

#3.4 Estimate Transcript Accessibility for the RNA Target at Single-Nucleotide Resolution
accessibility = RNAplfold_R(RNA_target,u.in = max(oligo_lengths)) %>%
  as_tibble() %>%
  mutate(end=1:l) %>%
  gather(length,accessibility, -end) %>%
  mutate(length=as.double(length))
target_annotation = left_join(target_annotation,accessibility,by=c('length','end'))


###########################################High-Frequency Polymorphisms analysis

#correcting end and start cord based on direction   
if (chr_coord['strand'] == 1) {
  target_annotation$chr_start = chr_coord['start'] + target_annotation$start - 1
  target_annotation$chr_end = chr_coord['start'] + target_annotation$end - 1
} else {
  target_annotation$chr_start = chr_coord['end'] - target_annotation$end + 1
  target_annotation$chr_end = chr_coord['end'] - target_annotation$start + 1
}

#keep unique names only and extract 
#information base on chr_start from target.
PM_freq = PMs %>%
  mutate(name = map(chr_start, function(X){
    filter(target_annotation,
           chr_start <= X,
           chr_end >= X) %>%
      select(name, chr_start_anno = chr_start)
  })) %>%
  unnest_legacy() %>%
  rename(chr_start_PM = chr_start) %>%
  group_by(name, chr_start_anno) %>%
  summarise(
    PM_tot_freq = 1 - prod(1 - PM_freq),
    PM_max_freq = max(PM_freq),
    PM_count= n()
  )

#combine data frames
target_annotation = left_join(
  target_annotation, PM_freq,
  by = c("name" = "name", "chr_start" = "chr_start_anno"))


##################################Match RNA Target Regions to the Mouse Ortholog

#get length
lm = width(RNA_target_mouse)

#make table of mouse information
MM_tab = lapply(oligo_lengths, function(i){
  tibble(st=1:(lm-i+1),w=i)}) %>%
  bind_rows()

#makes dnastringset object with mouse info.
RNAsitesMM = DNAStringSet(
  RNA_target_mouse[[1]],
  start=MM_tab$st,
  width=MM_tab$w)

#adds if conserved in mouse.
uni_tar$conserved_in_mmusculus = uni_tar$name %in% RNAsitesMM

################Occurrences of CGs Nucleobase Sequences, Duplex and Self-Folding

#bereken aantal "cg"
uni_tar$CGs = ( uni_tar$length -
                  nchar(gsub('CG','',uni_tar$name)) )/2
#maak de unique reverse complements van sequences in target regions
nucleobase_seq =unique(reverseComplement(target_regions))

#voeg ze toe aan uni_tar gekoppeld aan name
uni_tar$oligo_seq = as.character(nucleobase_seq[uni_tar$name])

#### deze nog toevoegen wanneer viennaRNA werkt
uni_tar$sec_energy = RNAselffold_R(uni_tar$oligo_seq)   #ram problem, make sure ram is free.  maybe use gc()
uni_tar$duplex_energy = RNAduplex_R(uni_tar$oligo_seq)  #ram problem, make sure ram is free.

#toxicity score acute neurotox score (desired >70)
uni_tar$tox_score = calculate_acute_neurotox(uni_tar$oligo_seq)



################################Selecting Nucleobase Sequences for Gapmer Design

#join tables
target_regions = left_join(
  target_annotation,
  uni_tar,
  by=c('name','length'))

#filter on needs >> this will all be options
target_region_select = filter(
  target_regions,
  gene_hits_pm ==1,
  gene_hits_1mm <= 50,
  accessibility > 1E-6,
  PM_tot_freq < 0.05,
  conserved_in_mmusculus,
  tox_score >= 70
  )

####################################################################################


library(cluster)
#get all unique starting positions in order
start_pos = sort(unique(target_region_select$start))

num_data_points <- length(start_pos)

#Set amount of clusters
if (num_data_points > 1) {
  K <- min(10, num_data_points - 1)
} else {
  # Handle the case when there are not enough data points
  K <- 1
}

print(K)

cluster_tab = tibble(
  start = start_pos,
  cluster = clara(x = start_pos, k = K,
                  metric = 'euclidean',
                  pamLike = T, samples=100)$clustering)
nucleobase_select = left_join(target_region_select,
                              cluster_tab, by='start') %>%
  group_by(cluster) %>%
  sample_n(1) %>%
  ungroup()



#collect the data, change the name for each gene tested
write_excel_csv(target_region_select,"SCN2A_selected.csv")
write_excel_csv(nucleobase_select,"SCN2A_final_10.csv")
