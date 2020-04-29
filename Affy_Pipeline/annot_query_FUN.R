## This function searches for the one gene symbol in the provided HGNC release file
annot_query_FUN <- function (gS, sym, eid, syn, psym = NA){
  if (length(grep(paste("^(", gS, ")$", sep = ""), sym)) > 0) {
    idx <- grep(paste("^(", gS, ")$", sep = ""), sym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")$", sep = ""), psym)) > 0) {
    idx <- grep(paste("^(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")\\|", sep = ""), psym)) > 0) {
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")\\|", sep = ""), psym)) > 0) {
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")$", sep = ""), psym)) > 0) {
    idx <- grep(paste("\\|(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (length(grep(paste("^(", gS, ")$", sep = ""), syn)) > 0) {
    idx <- grep(paste("^(", gS, ")$", sep = ""), syn)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else {
    idx <- grep(paste("^(", gS, ")\\|", sep = ""), syn)
    idx2 <- grep(paste("\\|(", gS, ")\\|", sep = ""), syn)
    idx3 <- grep(paste("\\|(", gS, ")$", sep = ""), syn)
    if (length(idx) > 0) {
      return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
    } else if (length(idx2) > 0) {
      return(c(paste(sym[idx2], collapse = ";"), paste(eid[idx2], collapse = ";")))
    } else {
      return(c(paste(sym[idx3], collapse = ";"), paste(eid[idx3], collapse = ";")))
    }
  }
}
