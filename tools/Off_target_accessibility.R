library(tidyverse)

acc_snippet <- function(
    begin_target,
    end_target,
    begin_snippet,
    end_snippet,
    full_snippet,
    total_length = 80
) {
  
  # This function extracts the snippet and target sequence information obtained from GGGenome 
  # for off-targets to create an 80-nt snippet with the off-target sequence as central as possible. 
  # The function returns a data frame containing the sequence of the 80-nt snippet, 
  # the snippet length (80 by default), the positions of the 80-nt snippet, and the target sequence positions.
  
  target_length <- end_target - begin_target + 1
  
  flank_total <- total_length - target_length
  flank_left  <- floor(flank_total / 2)
  flank_right <- ceiling(flank_total / 2)
  
  snippet_start <- begin_target - flank_left
  snippet_end   <- end_target   + flank_right
  
  if (snippet_start < begin_snippet) {
    shift <- begin_snippet - snippet_start
    snippet_start <- begin_snippet
    snippet_end   <- snippet_end + shift
  }
  
  if (snippet_end > end_snippet) {
    shift <- snippet_end - end_snippet
    snippet_end   <- end_snippet
    snippet_start <- snippet_start - shift
  }
  
  target_start_internal <- begin_target - snippet_start + 1
  target_end_internal   <- end_target   - snippet_start + 1
  
  start_i <- snippet_start - begin_snippet + 1
  end_i   <- snippet_end   - begin_snippet + 1
  
  snippet_seq <- substr(full_snippet, start_i, end_i)
  
  return(
    list(
      snippet_seq             = snippet_seq,
      snippet_start_external  = snippet_start,
      snippet_end_external    = snippet_end,
      snippet_length          = total_length,
      target_start_internal   = target_start_internal,
      target_end_internal     = target_end_internal
    )
  )
}
