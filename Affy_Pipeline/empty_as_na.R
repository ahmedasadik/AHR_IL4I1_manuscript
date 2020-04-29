## This function turns blank fields into NA
empty_as_na <- function (x) {
  if ("factor" %in% class(x)) 
    x <- as.character(x)
  ifelse(as.character(x) != "", x, NA)
}