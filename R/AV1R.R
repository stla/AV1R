#' Sum of squares
#' 
#' @export
#' @import dplyr
SOS <- function(y, group){
  return(data.frame(y=y,group=group) %>% mutate(mean=mean(y)) %>% group_by(group) %>% mutate(groupmean=mean(y)) %>% ungroup %>% do(summarise(., SSw=crossprod(y-groupmean), SSb=crossprod(groupmean-mean))) %>% as.list)
}

#' Simulates balanced design
#' 
#' @export
#' @importFrom magrittr "%>%"
SimAV1balanced <- function(I=3, J=4, mu=0, sigmab=2, sigmaw=1, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  return(data.frame(sapply(rnorm(I,mu,sigmab), function(mu) rnorm(J,mu,sigmaw))) %>%
  setNames(LETTERS[1:I]) %>% stack %>% setNames(c("y", "group")))
}

#' Samples Box & Tiao posterior distributions
#' 
#' @export
#' @import data.table
ranova <- function(y, group, nsims=50000){
  group <- factor(group)
  I <- length(levels(group))
  if(!all(I*as.numeric(table(group))==length(y))){ stop("balanced only!") }
  J <- length(y)/I
  SS <- SOS(y, group) #  sum of squares SSb SSw
  SSb <- SS$SSb; SSw <- SS$SSw
  omean <- mean(y)
  k <- tot <- 0
  z <- rnorm(nsims)
  ETAb <- ETAw <- rep(NA,nsims)
  while(k < nsims){
    tot <- tot +1
    etab <- SSb/rchisq(1, I-1)  
    etaw <- SSw/rchisq(1, I*(J-1))
    if(etab>etaw){
      k <- k+1
      ETAb[k] <- etab
      ETAw[k] <- etaw
    }
  }
  out <- data.table(mu=omean + sqrt(ETAb/I/J)*z, sigma2w=ETAw, sigma2b=(ETAb-ETAw)/J)
  reject <- 100*(tot-nsims)/tot
  cat("Rejection percentage: ", round(reject,1), "%")
  setattr(out, "reject", reject)
  return(out)
}