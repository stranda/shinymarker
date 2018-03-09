#
# take a group of directories and find all the fsa files.
# run the analyses and report the filenames 
#

find.files <- function(
    folders <- c("/home/astrand/GoogleDrive/data/sparrow/raw_fsa/Jan2014/"), #vector of folders
    xlim.size = c(75,250) #sizes to save for binary object.  No loci should have alleles that exceed....
    )
{
    everybody <- lapply(folders,function(fld)
                        {
                            files <- list.files(path=fld,"*.fsa")
                            do.call(rbind,
                                    mclapply(files,mc.cores=3,function(fn)
                                             {
                                                 contents=tryCatch({analyze.standard(read.fsa(files=fn,path=fld,sig.channel=1:4))},
                                                     warning=function(war){paste("war")},
                                                     error=function(err){NULL})
                                                 if (!is.null(contents)) {contents$data <- NULL; qual <- summary(contents$model)$r.squared} else
                                                     {qual <- 0}
                                                 data.frame(path=fld,
                                                            file=fn,
                                                            name=strsplit(fn,"\\.")[[1]][1],
                                                            qual=qual
                                                            )
                                                 
                                             })
                                    )
                        })
    do.call(rbind,everybody)
}
#all.data <- everybody
#save(file="all-data.rda",all.data)
    
