
#liz.500 <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)
analyze.standard <- function(fsadat,
                             actual.peaks=c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                             min.signal=500,
                             xlim=c(1200,7000),  # range of times
                             order=3
                             )
    {
        if (FALSE) { actual.peaks=c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500);min.signal=500;  xlim=c(1200,7000)}
        tmp <- fsadat
        ## first clip the sequence times to the raw xlimits
        tmp <- tmp[(tmp$time>=xlim[1])&(tmp$time<=xlim[2]),]
        ## make a copy of just the standard data 
        std <- tmp[tmp$chan=="standard",]
        ##reverse the direction of the sequence
        #drop out a lower threshold on standard
        std$peak <- ifelse(std$peak<min.signal,min.signal,std$peak)

        #now, try to find the peaks.  Using an approach from stackexchange
        peaks <- as.zoo(std$peak)
        peak.times <- sort(-1 * std$time[which(rollapply(peaks, 3, function(x) which.max(x)==2) )])
        if (length(peak.times)>length(actual.peaks))
            {
                peak.times <- peak.times[1:length(actual.peaks)]
            }
        if (length(peak.times)<length(actual.peaks))
            {
                actual.peaks <- sort(sort(actual.peaks,decreasing=TRUE)[1:length(peak.times)])
            }

        peak.times <- sort(abs(peak.times))

        std.curve <- lm(actual.peaks~poly(peak.times,as.numeric(order))) 
        tmp$size <- predict(std.curve,newdata=data.frame(peak.times=tmp$time))
        sse.predict <- sum((predict(std.curve)-actual.peaks)^2)
        list(data=tmp, model=std.curve,sse=sse.predict)
    }
