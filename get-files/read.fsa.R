read.fsa <- function(files = NULL, path = "./", sig.channel = 1:3, lad.channel = 105, pretrim = FALSE,
                     posttrim = ".fsa", thresh = -100, verbose = TRUE){

    require(seqinr)
  if(is.null(files))
    files <- list.files(path, pattern = "\\.fsa$", full.names = TRUE)
  else
    files <- paste(path, files, sep = "")

  res <- do.call(rbind, lapply(files, function(file) {
    if (verbose) message(file)
    abif <- read.abif(file)
    tag <- tag.trimmer(basename(file), pretrim, posttrim)
    
    lad.dat <- abif$Data[[paste('DATA.', lad.channel, sep='')]]
    
    res1 <- data.frame(tag = as.character(rep(tag, length(lad.dat))),
                       chan = as.character(rep("standard", length(lad.dat))),
                       time = as.numeric(1:length(lad.dat)),
                       peak = as.numeric(lad.dat))
    
    for (i in sig.channel) {
      chan.dat <- abif$Data[[paste('DATA.', i, sep='')]]
      res1 <- rbind(res1, data.frame(tag = as.character(rep(tag, length(chan.dat))),
                                     chan = as.character(rep(i, length(chan.dat))),
                                     time = as.numeric(1:length(chan.dat)),
                                     peak = as.numeric(chan.dat)))
    }
    res1
  }))
    
  if (thresh > -10) res <- subset(res, peak > thresh)
  return(res)
}

tag.trimmer <- function(x, pretrim = FALSE, posttrim = FALSE) {
  if(! is.na(pretrim)) {
    x <- sub(paste("^", pretrim, sep = ""), "", x)
  }
  if(! is.na(posttrim)){
    x <- sub(paste(posttrim, "$", sep = ""), "", x)
  }
  x
}
