#
# take a group of directories and find all the fsa files.
# run the analyses and report the filenames 
#

find.files <- function(
    folders= c("/home/astrand/GoogleDrive/data/sparrow/raw_fsa/Jan2014/"), #vector of folders
    xlim.size = c(75,500) #sizes to save for binary object.  No loci should have alleles that exceed....
)
{
    everybody <- lapply(folders,function(fld)
        {
            files <- list.files(path=fld,"*.fsa")
            do.call(rbind,
                    lapply(files,function(fn)
                        {
                            contents=tryCatch({analyze.standard(read.fsa(files=fn,path=fld,sig.channel=1:4))},
                                warning=function(war){paste("war")},
                                error=function(err){NULL})
                            if (!is.null(contents)) {contents$data <- NULL; } 
                            
                            qual <- contents$sse
                            print(fld)
                            print(fn)
                            print(strsplit(fn,"\\.")[[1]][1])
                            print(qual)
                            tryCatch({data.frame(path=fld,
                                       file=fn,
                                       name=strsplit(fn,"\\.")[[1]][1],
                                       qual=qual
                                       )},
                                     error=function(error){data.frame(path=fld,file=fn,name=strsplit(fn,"\\.")[[1]][1],qual=Inf)})
                            
                        })
                    )
        })
    all.files <- do.call(rbind,everybody)
    save(file="file.status.rda",all.files)
}


#all.data <- everybody
#save(file="all-data.rda",all.data)
    
