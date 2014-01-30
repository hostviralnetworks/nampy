# Author: ghannum
###############################################################################
library(actuar)
library(MASS)


fitPareto <- function(x, attempt=0)
{
    val <- c()
    
    ft <- function(attempt)
    {
            if (attempt==0) val <- suppressWarnings(fitdistr(x,densfun=dpareto,start=list(shape=300,scale=10))$estimate)
            else if (attempt==1) val <- suppressWarnings(fitdistr(x,densfun=dpareto,start=list(shape=3000,scale=1000))$estimate)
            else if (attempt==2) val <- suppressWarnings(fitdistr(x,densfun=dpareto,start=list(shape=1000,scale=1))$estimate)
            else return(NULL)
    }
    
    val <- tryCatch(ft(attempt),error=function(e){NULL})
            
    if (is.null(val) && attempt<2)    return( fitPareto(x, attempt+1) )
    
    if (is.null(val)) warning("Pareto fit failed")
    
    return(val)
}


getParetoFit <- function(d, side="upper", tailSize=50, na.rm=TRUE)
{
    d <- na.omit(d)
    
    d <- sort(d)
    N <- length(d)
    
    pfit1 <- NULL
    pfit2 <- NULL
    
    if (side=="upper" || side=="two-sided")
    {
            tail <- d[(N-tailSize):N]
            
            t1 <- (d[N-tailSize]+d[N-tailSize-1])/2
            
            tail <- tail - t1
            
            tol <- 1e-10
            index = 2;
            while (abs(tail[index]-tail[index-1])<tol)
            {
                    tail[index] <- tol;
                    index <- index + 1
                    
                    if (index>=length(tail))
                    {
                            warning("All-same tail1")
                            return(NULL)
                    }
            }

            pfit1 <- fitPareto(tail)
    }
    
    if (side=="lower" || side=="two-sided")
    {
            tail <- d[tailSize:1]
            
            t2 <- (d[tailSize]+d[tailSize+1])/2
            
            tail <- abs(tail - t2)
            
            tol <- 1e-10
            index = 2;
            while (abs(tail[index]-tail[index-1])<tol)
            {
                    tail[index] <- tol;
                    index <- index + 1
                    
                    if (index>=length(tail))
                    {
                            print(tail)
                            warning("All-same tail2")
                            return(NULL)
                    }
            }
            
            pfit2 <- fitPareto(tail)
    }
    
    
    ebounds <- c(d[11],d[N-11])
    cdf <- ecdf(d)
    
    if(is.null(pfit1)) pfit1 <- c(NA,NA)
    if(is.null(pfit2)) pfit2 <- c(NA,NA)
    
    if (side=="upper") return(list(cdf=cdf,shape1=pfit1[1],scale1=pfit1[2],t1=t1,shape2=NULL,scale2=NULL,t2=NULL,tailFrac=tailSize/length(d),side=side,ebounds=ebounds,N=length(d)))
    else if (side=="lower") return(list(cdf=cdf,shape1=NULL,scale1=NULL,t1=NULL,shape2=pfit2[1],scale2=pfit2[2],t2=t2,tailFrac=tailSize/length(d),side=side,ebounds=ebounds,N=length(d)))
    else return(list(cdf=cdf,shape1=pfit1[1],scale1=pfit1[2],t1=t1,shape2=pfit2[1],scale2=pfit2[2],t2=t2,tailFrac=tailSize/length(d),side=side,ebounds=ebounds,N=length(d)))
}


paretoExtrapolate <- function(x, fit)
{
    if (length(x)>1) return( sapply(x, paretoExtrapolate, fit=fit) )
    
    if (is.na(x)) return(NA)
    
    if (fit$side=="upper")
    {
            if (x > fit$ebounds[2] && is.na(fit$shape1)) warning("Failed fit, using empirical p-value")
            if (x > fit$ebounds[2] && !is.na(fit$shape1))    return(ppareto(x-fit$t1, fit$shape1, fit$scale1, lower.tail=FALSE)*fit$tailFrac)
            else return(max(1-fit$cdf(x),.5/fit$N))
    }else if (fit$side=="lower")
    {
            if (x < fit$ebounds[1] &&is.na(fit$shape2)) warning("Failed fit, using empirical p-value")
            if (x < fit$ebounds[1] && !is.na(fit$shape2))    return(ppareto(fit$t2-x, fit$shape2, fit$scale2, lower.tail=FALSE)*fit$tailFrac)
            else return(max(fit$cdf(x),.5/fit$N))
    }else
    {
            if ((x < fit$ebounds[1] && is.na(fit$shape2)) || (x > fit$ebounds[2] && is.na(fit$shape1))) warning("Failed fit, using empirical p-value")
            if (x < fit$ebounds[1] && !is.na(fit$shape2))    return(2 * ppareto(fit$t2-x, fit$shape2, fit$scale2, lower.tail=FALSE)*fit$tailFrac)
            else if (x > fit$ebounds[2] && !is.na(fit$shape1)) return(2 * ppareto(x-fit$t1, fit$shape1, fit$scale1, lower.tail=FALSE)*fit$tailFrac)
            else
            {
                    p <- fit$cdf(x)
                    return(max(2 * min(p,1-p),1/fit$N))
            }
    }
}