## This function removes Vs from the file names of transcription factor target files
Vs_FUN <- function (x){
  if (length(grep("^V\\_", x)) > 0) {
    f_n <- gsub("V\\_", "", x) %>% gsub("_.*", "", .) %>% 
      gsub("\\..*", "", .)
  } else {
    f_n <- x %>% gsub("_.*", "", .) %>% gsub("\\..*", "", .)
  }
}
