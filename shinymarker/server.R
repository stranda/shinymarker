source("setup.R")
#### Define server logic required to summarize and view the selected dataset
updateflag <<- 0
hoseProjectflag <<- 0

shinyServer(function(input, output) {
    values <- reactive({
        names <- input$inds

#        comp <- read.csv("comp.csv")
#        alleles <- list(mme1=as.numeric(unique(c(comp$MME1.a,comp$MME1.b))),
#                        mme2=as.numeric(unique(c(comp$MME2.a,comp$MME2.b))),
#                        mme3=as.numeric(unique(c(comp$MME3.a,comp$MME3.b))),
#                        mme7=as.numeric(unique(c(comp$MME7.a,comp$MME7.b))),
#                        mme8=as.numeric(unique(c(comp$MME8.a,comp$MME8.b))),
#                        mme12=as.numeric(unique(c(comp$MME12.a,comp$MME12.b))))
             
                                        #        names <- strsplit(gsub("\ ","",input$inds),",")[[1]]
        ret <- lapply(names,function(nm)
                      {
                          fn <- paste0(nm,".fsa")
                          if (file.exists(paste0(basepath,input$project,"/",fn)))
                              {
                                  tmp <- analyze.standard(read.fsa(files=fn,path=paste0(basepath,input$project,"/"),
                                                                   sig.channel=1:4),
                                                          xlim=c(input$min.time,9000),
                                                          min.signal=input$min.signal,
                                                          order=input$order)
                                  ret <- list(model= tmp$model,
                                              data= tmp$data,
                                              sse=tmp$sse,
                                              path=path,
                                              file= fn,
                                              name=nm,
                                              alleles=NULL)
                                  ret
                              } else {NULL}
                      })
        
    })
    
    observe({
        projpath <- paste0(basepath,isolate(input$project))
#        if (!file.exists(projpath)) {dir.create(projpath)}
        inFile <- input$filename

        if ((!is.null(inFile)))
            {
                if (file.exists(inFile$datapath))
                    unzip(inFile$datapath, junkpaths=TRUE, exdir=projpath)
            }
        if (input$hoseProject>hoseProjectflag)
            {
                hoseProjectflag <<- input$hoseProject
                unlink(inFile$datapath,force=T)
                unlink(projpath,recursive=T,force=T)
                
            }
        if (input$updateIndividuals>updateflag)
            {
                updateflag <<- input$updateIndividuals
                files <- list.files(path=paste0(basepath,input$project),pattern=".fsa")
                names <- c(sapply(strsplit(files,"\\."),function(x){x[[1]][1]}))
                 output$chooseind <- renderUI({
                   selectInput('inds',"Individuals",choices=names,selected=names[1],multiple=T)
                })
            }
    })


    
    output$size <- renderText(
        {
            paste("Allele size at mouse cursor",round(as.numeric(input$plot_hover$x),2))
        })
    
    output$test <- renderTable(
        {
            dat <- values()
            dat[[1]]
        })
    
    output$aligned.dye <- renderPlot(
        {
            dat <- values()
            n <- length(dat)
            dye <- as.numeric(as.character(input$dye))
            
            grid <- seq(0,500,as.numeric(as.character(input$grid)))
            par(mfrow=c(n+1,1))
            for (i in 1:length(dat))
                {
                    x <- dat[[i]]
                    xz <- input$xzoom
                    center <- input$pointer
                    xlim <- c(center-(xz/2),center+(xz/2))
                    ylim <- max(x$data$peak)*input$yrange
                    if (i==1)
                        {
                            par(mar=c(0,4,2,1))
                            plot(peak~size,type="n",dat=x$data,
                                 xlim=xlim,ylim=ylim,axes=F,xlab="",ylab="Signal",cex.lab=1.3)
                            axis(3,cex.axis=1.5)
                        } else {
                            if (i<length(dat))
                                {
                                    par(mar=c(0,4,0,1))
                                    plot(peak~size,type="n",col=1,dat=x$data,
                                         xlim=xlim,ylim=ylim,axes=F,xlab="",ylab="Signal",cex.lab=1.3)
                                } else {
                                    par(mar=c(1,4,0,1))
                                    plot(peak~size,type="n",dat=x$data,
                                         xlim=xlim,ylim=ylim,axes=F,xlab="Size",ylab="Signal",cex.lab=1.3)
                                    axis(1,cex.axis=1.5)
                                }
                        }
                    box()
                    axis(2)
                    
                    abline(v=grid,col="lightgray")
                    
                    for (j in dye)
                        {
                                        #                            if (input$show.alleles) {
                                        #                            acrd <- unlist(x$alleles[c(unlist(list(c(4,5),c(1,6),c(2),c(3))[[j]]))])
                                        #                            abline(v=acrd,col=j+1,lwd=2)
                                        #                        }
                            if (j<4) col=j+1 else col='darksalmon'
                            points(peak~size,type="l",col="white",lwd=2,dat=x$data[x$data$chan==j,])
                            points(peak~size,type="l",col=col,dat=x$data[x$data$chan==j,])
                        }
                    
                    if (input$standard) {points(peak~size, type="l",
                                                col="magenta", dat=x$data[x$data$chan=="standard",])}
                    
                    text(x=xlim[1],y=0.8*ylim[2],paste("id",x$name),cex=1.4,pos=4)
                    abline(v=center,lwd=1.5,col="red")
                    text(x=center,y=0.9*ylim[2],pos=2,center,cex=2)
                                        #                            rsq <- summary(x$model)$r.squared
                    sse <- x$sse
                                        #                            text(x=xlim[1],y=0.7*ylim[2], round(rsq,3),pos=4)
                    
                    text(x=xlim[2],y=0.9*ylim[2], paste("Qual (lower better)",round(sse,1)),col=ifelse(sse<30,"green",ifelse(sse>80,"red","orange")),pos=2,cex=1.5)
                }
            
        },height=function(){length(input$inds)*150})
    
    output$sameplot <- renderPlot(
        {
            dat <- values()
            n <- length(dat)
            dye <- as.numeric(as.character(input$dye))
            
            grid <- seq(0,500,as.numeric(as.character(input$grid)))
            xz <- input$xzoom
            center <- input$pointer
            xlim <- c(center-(xz/2),center+(xz/2))
            
            ylim <- max(sapply(dat,function(x){max(x$data$peak)})) * input$yrange
            
            plot(peak~size,type="n",dat=dat[[1]]$data,
                 xlim=xlim,ylim=ylim,axes=T,xlab="Size (bp)",ylab="Signal",cex.lab=1.3)
            
            abline(v=grid,col="lightgray")
            abline(v=center,lwd=1.5,col="red")
            text(x=center,y=0.9*ylim[2],pos=2,center,cex=2)
            txt <- NULL
            for (i in 1:length(dat))
                {
                    x <- dat[[i]]
                    txt <- c(txt,x$name)
                    for (j in dye)
                        {
                                        #                            if (input$show.alleles) {
                                        #                            acrd <- unlist(x$alleles[c(unlist(list(c(4,5),c(1,6),c(2),c(3))[[j]]))])
                                        #                            abline(v=acrd,col=j+1,lwd=2)
                                        #                        }
                            points(peak~size,type="l",col="white",lwd=3,lty=j,dat=x$data[x$data$chan==j,])
                            points(peak~size,type="l",col=i+1,lty=j,dat=x$data[x$data$chan==j,])
                        }
                    if (input$standard) {points(peak~size, type="l", col="magenta", dat=x$data[x$data$chan=="standard",])}
                    
                }
            legend(x=xlim[1],y=ylim[2]*0.95,legend=txt,col=1+(1:length(dat)),lty=1)
            
        })
    
})
