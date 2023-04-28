entropy_O1 <- function(SEQ)
{

  # this function computes for any sequence SEQ, its Shannon entropy, and the corrected conditional entropies at orders 1 to 5 following Porta et al. (1998)
  # the output is: 
          # O_corr: the vector of the Shannon entropy followed by the corrected conditional entropies for order 1 to 5
  
  # max4 and max5 are the maximum numbers of different sites visited in SEQ for the conditional entropies of order 4 and 5 to be computed, respectively
  # memory size problems can arise when computing the highest-order conditional transition matrices, for sequences with a high number of different sites used.
       # in this case, an error message is given, you can try to increase the amount of memory used by R using memory.limit()
       # otherwise, max4 and max5 should be set to smaller values
  # if the number of different sites visited is superior to max3, the computation is stopped and returns a series of NA
  
  # In any case, the relevant estimation of the minimum conditional entropy is warranted only if a minimum is reached in O_corr
  # that is, that the conditional entropy does not significantly decrease until the last value (that is, O5, or O4 or O3 if a very high number of different sites are visited)

  
  sites <- unique(SEQ)
  base <- length(sites)

  #-------------------------------
  ## order 0
  #-------------------------------

  prop <- sapply(sites,function(x) length(SEQ[SEQ==x])/length(SEQ))

  O0 <- -sum(sapply(1:base,function(x) prop[x]*log(prop[x],base)))

  #------------------------------
  ## order 1
  #------------------------------

  p1 <- matrix(nrow = base, ncol = base, 0)      # counts the transitions from one site to another
  for (t in 1:(length(SEQ) - 1)) p1[which(sites==SEQ[t]), which(sites==SEQ[t + 1])] <- p1[which(sites==SEQ[t]), which(sites==SEQ[t + 1])] + 1
  
  prop1 <- matrix(nrow = base, ncol = base, 0)   # conditional probabilities
  for (i in 1:base) 
    { if(sum(p1[i,])>0){prop1[i, ] <- p1[i, ] / sum(p1[i, ])   }}

  O1 <- 0
  for(i in 1:base)
  {
    for(j in 1:base)
    {
      if(prop1[i,j]>0)
      {O1 <- O1 - prop[i]*prop1[i,j]*log(prop1[i,j],base)}
    }
  }

  O1_corr <- O1 + (length(which(p1==1))/(length(SEQ)-1))*O0

  rm(prop1)
  
  return(min(O0, O1_corr))

}
