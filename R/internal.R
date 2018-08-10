globalVariables(c("log_linear_match","log_linear_match_c","log_linear_mismatch","log_linear_mismatch_c"))

getnewindexes_all <- function(newterm,colnames){
  indexes <- sort(unlist(strsplit(newterm,split=":",fixed=TRUE)))
  toadd <- setdiff(colnames,indexes)
  newterm_big <- character(length(toadd))
  for(i in 1:length(toadd)){
    newterm_big[i] <- paste(sort(c(indexes,toadd[i])),c(rep(":",length(indexes)),""),sep="",collapse="")
  }
  newterm_big
}

###  download mortality matrix

find_pattern <- function(count_dframe, vec){
  indexes <- which(vec %in% c(0,1))
 
  if(length(indexes)==0){
    rows = 1:nrow(count_dframe)
  }else{
    subvec <- vec[indexes]
    ifelse(length(indexes)==1, rows <- which(sapply(count_dframe[,indexes],function(d) all(d == subvec))),
           rows <- which(apply(count_dframe[,indexes],1,function(d) all(d == subvec))))
    }
  rows
}

create_01mat = function(L){
  ncombin <- 2^L
  
  mat = data.frame(v=rep.int(c(0,1),rep.int(ncombin/2,2)))
  for(coln in 2:L){
    mat <- cbind(mat, v=rep.int(rep.int(c(0,1),rep.int(ncombin/(2^coln),2)),2^(coln-1)))
    
  }
  mat
}

##  EM algorithm to fill in matrix

imputemissing <- function(count_dframe,zero_one_dframe,tol=10^-4){
  nc <- ncol(zero_one_dframe)
  
  flevs = 1:nrow(zero_one_dframe)
  indices = apply(count_dframe[,1:nc],1,function(d) find_pattern(zero_one_dframe,d))
  indexes = factor(unlist(indices),levels=flevs)
  oldprobs <- rep.int(1/nrow(zero_one_dframe),nrow(zero_one_dframe))
  newcounts = unlist(mapply(function(counts,i) oldprobs[i]/sum(oldprobs[i])*counts,
                            count_dframe[,nc+1],indices))
  
  Ecounts = tapply(newcounts,indexes,sum,default = 0)
  newprobs <- Ecounts/sum(Ecounts)
  
  while(max(abs(newprobs - oldprobs))>tol){
    oldprobs <- newprobs
    
    newcounts = unlist(mapply(function(counts,i) oldprobs[i]/sum(oldprobs[i])*counts,count_dframe[,nc+1],indices))
    
    Ecounts = tapply(newcounts,indexes,sum,default = 0)
    
    newprobs <- Ecounts/sum(Ecounts)
  }
  return(list(counts=Ecounts,probs=newprobs))
}




###  reassign probabilities for matrix which has 2's.

reassign_probs <- function(comparemat, zero_one_dframe, probabilities){
  nc <- ncol(comparemat)
  zodf = zero_one_dframe[,1:nc]
  counts = zero_one_dframe$counts
  
  indices = apply(comparemat,1,function(v) find_pattern(zodf,v))
  sapply(indices, function(i) sum(counts[i] * probabilities[i])) / sapply(indices, function(i) sum(counts[i]))
}

### mortality ###


EM_match_modelsearch <- function(count_dframe,m,u,fixedcol=c(2:9),p_init=0.5,tol=10^-5,maxit=10000,allterms=allterms){
  N <- sum(count_dframe$counts)
  L <- ncol(count_dframe)-1
  
  stuff <- EM_match_independence_v3(count_dframe,m,u,p_init=0.5,tol=10^-5,fixedcol=fixedcol)  ##  get starting values.
  p_new <- stuff$p
  p_old <- p_new
  gamma_current <- stuff$probs
  counter=1
  probs_match_new <- count_dframe$counts*gamma_current/sum(count_dframe$counts*gamma_current)
  probs_match_old <- probs_match_new
  probs_mismatch_new <- count_dframe$counts*(1-gamma_current)/sum(count_dframe$counts*(1-gamma_current))
  probs_mismatch_old <- probs_mismatch_new
  
  Ecounts_match = count_dframe
  Ecounts_mismatch = count_dframe
  
  while((counter<=maxit & max(abs(c(probs_match_new-probs_match_old, probs_mismatch_new-probs_mismatch_old,p_new-p_old))) > tol) |counter<=2){
    
    probs_match_old <-  probs_match_new
    probs_mismatch_old <-  probs_mismatch_new
    p_old <- p_new
    
    gamma_current <- p_old*probs_match_old/(p_old*probs_match_old  + (1-p_old)*probs_mismatch_old)
    gamma_current[is.na(gamma_current)] <- 0
    
    p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)
    
    Ecounts_match$counts = count_dframe$counts*gamma_current
    Ecounts_mismatch$counts = count_dframe$counts*(1-gamma_current)
    
    ##  only do the following on first iteration:
    if(counter==1){
      
      newterm = 1
      logN <- log(sum(count_dframe$counts))
      
      termsnotadded = allterms
      form_m = 'counts ~ .'
      log_linear_match <- glm(form_m,quasipoisson , Ecounts_match)
      BIC_0 = log_linear_match$deviance
      minBIC = BIC_0 / 2
      
      while(BIC_0 > minBIC){
        form_m = paste(form_m,newterm,sep='+')
        
        log_linear_match <- glm(form_m ,quasipoisson , Ecounts_match)
        
        dfres_0 = log_linear_match$df.residual
        BIC_0 = log_linear_match$deviance
        
        llmcs = lapply(termsnotadded,function(t){
          formc = paste(form_m, t, sep = '+')
          glm(formc,quasipoisson , Ecounts_match)
        })
        
        ds = sapply(llmcs,deviance)
        dfrs = sapply(llmcs, df.residual)
        
        BICvec <- ds + logN * (dfres_0 - dfrs)
        
        newterm = termsnotadded[which.min(BICvec)]
        
        termsnotadded = setdiff(termsnotadded,newterm)
        termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,colnames(count_dframe)[1:L]))
        
        minBIC = min(BICvec) 
      }
      
      newterm = 1
      
      termsnotadded = allterms
      form_u = 'counts ~ .'
      log_linear_mismatch <- glm(form_u,quasipoisson , Ecounts_mismatch)
      BIC_0 = log_linear_mismatch$deviance
      minBIC = BIC_0 / 2
      
      while(BIC_0 > minBIC){
        form_u = paste(form_u,newterm,sep='+')
        
        log_linear_mismatch <- glm(form_u ,quasipoisson , Ecounts_mismatch)
        
        dfres_0 = log_linear_mismatch$df.residual
        BIC_0 = log_linear_mismatch$deviance
        
        llmcs = lapply(termsnotadded,function(t){
          formc = paste(form_u, t, sep = '+')
          glm(formc,quasipoisson , Ecounts_mismatch)
        })
        
        ds = sapply(llmcs,deviance)
        dfrs = sapply(llmcs, df.residual)
        
        BICvec <- ds + logN * (dfres_0 - dfrs)
        
        newterm = termsnotadded[which.min(BICvec)]
        
        termsnotadded = setdiff(termsnotadded,newterm)
        termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,colnames(count_dframe)[1:L]))
        
        minBIC = min(BICvec) 
      }
    }
    if(counter > 1){
      log_linear_match <- glm(form_m,quasipoisson , Ecounts_match)
      log_linear_mismatch <- glm(form_u,quasipoisson , Ecounts_mismatch)
    }
    
    probs_match_new <- log_linear_match$fitted/sum(log_linear_match$fitted)
    probs_mismatch_new <-  log_linear_mismatch$fitted/sum(log_linear_mismatch$fitted)
    counter=counter + 1
  }
  errs = c(mean((Ecounts_match$counts - exp(predict(log_linear_match)))^2/exp(predict(log_linear_match))),
           mean((Ecounts_mismatch$counts - exp(predict(log_linear_mismatch)))^2/exp(predict(log_linear_mismatch))))
  nhats = c(sum(exp(predict(log_linear_match))), sum(exp(predict(log_linear_mismatch))))
  npar = c(length(coef(log_linear_match))-1, length(coef(log_linear_match))-1)
  alphas = c(p_new,1-p_new)
  
  MRC = sum(nhats*log(errs)) + sum(nhats*(nhats + npar)/(nhats-npar-2)) - 2*sum(nhats*log(alphas))
  
  return(list(p=p_new, probs = gamma_current, model_match = log_linear_match, model_mismatch=log_linear_mismatch,MRC=MRC))
}


EM_match_independence_v3 <- function(count_dframe,m,u,fixedcol=c(2:9),p_init=0.5,tol=10^-5){
  counts = count_dframe$counts
  N <- sum(counts)
  L <- ncol(count_dframe)-1
  u_col = (1:L)
  if(length(fixedcol) > 0) u_col = u_col[-fixedcol]
  
  gamma_current <- numeric(nrow(count_dframe))
  count_pattern_match <- numeric(nrow(count_dframe))
  count_pattern_nomatch <- numeric(nrow(count_dframe))
  
  
  #  start off the algorithm
  p_old <- p_init
  m_old <- m
  u_old <- u
  
  gamma_current <- apply(count_dframe[,1:L],1,function(v){
    p_init * prod(m_old^v) * prod((1-m_old)^(1-v)) / (p_init * prod(m_old^v) * prod((1-m_old)^(1-v)) + (1-p_init) * prod(u_old^v) * prod((1-u_old)^(1-v)))
  })
  
  p_new <- sum(gamma_current*counts)/N
  
  #create data frames for counts conditionally for both matches and non matches
  probs_match = counts*gamma_current/sum(counts*gamma_current)
  probs_mismatch = counts*(1-gamma_current)/sum(counts*(1-gamma_current))
  
  #Maximization
  m_new = colSums(count_dframe[,1:L] * probs_match)
  u_new = u 
  u_new[u_col] = colSums(count_dframe[,u_col] * probs_mismatch)
  
  while(max(abs(c(m_new-m_old, u_new-u_old,p_new-p_old))) > tol){
    p_old <- p_new
    m_old <- m_new 
    u_old <- u_new
    
    gamma_current <- apply(count_dframe[,1:L],1,function(v){
      p_old * prod(m_old^v) * prod((1-m_old)^(1-v)) / (p_old * prod(m_old^v) * prod((1-m_old)^(1-v)) + (1-p_old) * prod(u_old^v) * prod((1-u_old)^(1-v)))
    })
    
    p_new <- sum(gamma_current*counts)/N
    
    #Maximization
    m_new = colSums(count_dframe[,1:L] * probs_match)
    u_new = u 
    u_new[u_col] = colSums(count_dframe[,u_col] * probs_mismatch)
    
    probs_match = counts*gamma_current/sum(counts*gamma_current)
    probs_mismatch = counts*(1-gamma_current)/sum(counts*(1-gamma_current))
  }
  
  log_linear_match = glm(I(counts * gamma_current) ~ ., count_dframe, family = quasipoisson())
  log_linear_mismatch = glm(I(counts * (1-gamma_current)) ~ ., count_dframe, family = quasipoisson())
  
  errs = c(mean((counts * gamma_current - exp(predict(log_linear_match)))^2/exp(predict(log_linear_match))),
           mean((counts * (1-gamma_current) - exp(predict(log_linear_mismatch)))^2/exp(predict(log_linear_mismatch))))
  nhats = c(sum(exp(predict(log_linear_match))), sum(exp(predict(log_linear_mismatch))))
  npar = c(length(coef(log_linear_match))-1, length(coef(log_linear_match))-1)
  alphas = c(p_new,1-p_new)
  
  MRC = sum(nhats*log(errs)) + sum(nhats*(nhats + npar)/(nhats-npar-2)) - 2*sum(nhats*log(alphas))
  
  return(list(p=p_new, probs = gamma_current, model_match = log_linear_match, model_mismatch=log_linear_mismatch,MRC=MRC))
}


EM_match_modelsearch_iu <- function(count_dframe,m,u,fixedcol=c(2:9),p_init=0.5,tol=10^-5,maxit=10000,allterms=allterms){
  N <- sum(count_dframe$counts)
  L <- ncol(count_dframe)-1
  
  stuff <- EM_match_independence_v3(count_dframe,m,u,p_init=0.5,tol=10^-5,fixedcol=fixedcol)         ##  get starting values.
  p_new <- stuff$p
  p_old <- p_new
  gamma_current <- stuff$probs
  counter=1
  probs_match_new <- count_dframe$counts*gamma_current/sum(count_dframe$counts*gamma_current)
  probs_match_old <- probs_match_new
  probs_mismatch_new <- count_dframe$counts*(1-gamma_current)/sum(count_dframe$counts*(1-gamma_current))
  probs_mismatch_old <- probs_mismatch_new
  
  Ecounts_match = count_dframe
  Ecounts_mismatch = count_dframe
  
  while((counter<=maxit & max(abs(c(probs_match_new-probs_match_old, probs_mismatch_new-probs_mismatch_old,p_new-p_old))) > tol) |counter<=2){
    
    probs_match_old <-  probs_match_new
    probs_mismatch_old <-  probs_mismatch_new
    p_old <- p_new
    
    gamma_current <- p_old*probs_match_old/(p_old*probs_match_old  + (1-p_old)*probs_mismatch_old)
    gamma_current[is.na(gamma_current)] <- 0
    
    p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)
    
    Ecounts_match$counts = count_dframe$counts*gamma_current
    Ecounts_mismatch$counts = count_dframe$counts*(1-gamma_current)
    
    ##  only do the following on first iteration:
    if(counter==1){
      
      newterm = 1
      logN <- log(N)
      
      termsnotadded = allterms
      form = 'counts ~ .'
      log_linear_match <- glm(form,quasipoisson , Ecounts_match)
      BIC_0 = log_linear_match$deviance
      minBIC = BIC_0 / 2
      while(BIC_0 > minBIC){
        form = paste(form,newterm,sep='+')
        
        log_linear_match <- glm(form,quasipoisson , Ecounts_match)
        dfres_0 = log_linear_match$df.residual
        BIC_0 = log_linear_match$deviance
        llmcs = lapply(termsnotadded,function(t){
          formc = paste(form, t, sep = '+')
          glm(formc,quasipoisson , Ecounts_match)
        })
        
        ds = sapply(llmcs,deviance)
        dfrs = sapply(llmcs, df.residual)
        
        BICvec <- ds + logN * (dfres_0 - dfrs)
        
        newterm = termsnotadded[which.min(BICvec)]
        
        termsnotadded = setdiff(termsnotadded,newterm)
        termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,colnames(count_dframe)[1:L]))
        
        minBIC = min(BICvec) 
      }
      log_linear_mismatch <- glm(counts ~ .,quasipoisson , Ecounts_mismatch)
    }
    if(counter > 1){
      log_linear_match <- glm(form,quasipoisson , Ecounts_match)
      log_linear_mismatch <- glm(counts ~ .,quasipoisson , Ecounts_mismatch)
    }
    
    probs_match_new <- log_linear_match$fitted/sum(log_linear_match$fitted)
    probs_mismatch_new <-  log_linear_mismatch$fitted/sum(log_linear_mismatch$fitted)
    counter=counter + 1
  }
  errs = c(mean((Ecounts_match$counts - exp(predict(log_linear_match)))^2/exp(predict(log_linear_match))),
           mean((Ecounts_mismatch$counts - exp(predict(log_linear_mismatch)))^2/exp(predict(log_linear_mismatch))))
  nhats = c(sum(exp(predict(log_linear_match))), sum(exp(predict(log_linear_mismatch))))
  npar = c(length(coef(log_linear_match))-1, length(coef(log_linear_match))-1)
  alphas = c(p_new,1-p_new)
  
  MRC = sum(nhats*log(errs)) + sum(nhats*(nhats + npar)/(nhats-npar-2)) - 2*sum(nhats*log(alphas))
  
  return(list(p=p_new, probs = gamma_current, model_match = log_linear_match, model_mismatch=log_linear_mismatch,MRC=MRC))
}


gen_data <- function(means,correl=.1,nsample=500000, missingprobs){
  nvar <- length(means)
  ##  use a multivariate normal model:
  covmat <- matrix(correl,nrow=nvar,ncol=nvar)
  diag(covmat) <- 1
  stuff <- eigen(covmat)
  sqr_sigma <- stuff$vectors%*%diag(sqrt(stuff$values))%*%t(stuff$vectors)
  out_norm <- means + sqr_sigma%*%matrix(rnorm(nvar*nsample),nrow=nvar,ncol=nsample)
  out_match <- out_norm > 0
  mode(out_match) <- "numeric"
  for(j in 1:nvar){
    missingindex <- (1:nsample)[as.logical(rbinom(n=nsample,size=1,prob=missingprobs[j]))]
    if(length(missingindex)>0 ) out_match[j,missingindex] <- "2"
  }
  therows <- out_match[1,]
  for(i in 2:nvar) therows <- paste(therows,out_match[i,],sep="_")
  stuff <- table(therows)
  temp_dframe <- cbind(matrix(0,nrow=length(stuff),ncol=nvar),as.numeric(stuff))
  for(i in 1:nrow(temp_dframe)){
    temp_dframe[i,1:nvar] <- as.numeric(unlist(strsplit(names(stuff[i]),split="_")))
  }
  return(temp_dframe)

}

#' Function to simulate agreement patterns for record pairs from a mixture model
#' @param cor_match correlation in 0/1 agreement fields for record pairs constituting the same individual (matches).
#' @param cor_mismatch correlation in 0/1 agreement fields for record pairs not constituting the same individual (true mismatches).
#' @param nsample number of record pairs to simulate.
#' @param pi_match the probability both records from a randomly selected record pair correspond to the same individual (a true match).
#' @param m_probs marginal probabilities of agreement for each field for record pairs constituting true matches.
#' @param u_probs marginal probabilities of agreement for each field for record pairs constituting true mismatches.
#' @param missingprobs probabilities that each field is missing on at least one record pair.
#' @return A matrix of the simulated agreement patterns, with final column equal to the count of each pattern.
#' @keywords internal
#' @export
#' @importFrom stats qnorm rbinom rnorm
#' @importFrom utils flush.console
#' @examples
#' m_probs <- rep(0.8,6)
#' u_probs <- rep(0.2,6)
#' means_match <- -1*qnorm(1-m_probs)
#' means_mismatch <- -1*qnorm(1-u_probs)
#' missingprobs <- rep(.2,6)
#' thedata <- do_sim(cor_match=0.2,cor_mismatch=0,nsample=10^4,
#' pi_match=.5,m_probs=rep(0.8,5),u_probs=rep(0.2,5),missingprobs=rep(0.4,5))
#' thedata
do_sim <- function(cor_match=.1,cor_mismatch=0,nsample=10^3,pi_match=.01,m_probs,u_probs,missingprobs){
  means_match <- -1*qnorm(1-m_probs)
  means_mismatch <- -1*qnorm(1-u_probs)
    data_match <- gen_data(means=means_match,correl=cor_match,nsample=nsample*pi_match,missingprobs=missingprobs)
  data_mismatch  <- gen_data(means=means_mismatch,cor_mismatch,nsample=nsample*(1-pi_match),missingprobs=missingprobs)   ###  assume no corelation in matching status for patients not matching
  names_match <- data_match[,1]
  names_mismatch <- data_mismatch[,1]
  if(ncol(data_match)>=3)
  for(i in 2:(ncol(data_match)-1)){
    names_match  <- paste(names_match,data_match[,i],sep="_")
    names_mismatch <- paste(names_mismatch,data_mismatch[,i],sep="_")
  }
  data_match <- cbind(names_match,data_match)
  colnames(data_match)[ncol(data_match)] <- "count_match"
  data_mismatch <- cbind(names_mismatch,data_mismatch)
  colnames(data_mismatch)[ncol(data_match)] <- "count_mismatch"
  alldata <- merge(data_match,data_mismatch,by=1,all.x=TRUE,all.y=TRUE,as.is=TRUE)
  alldata$totalcount <- apply(cbind(as.numeric(as.character(alldata$count_match)),as.numeric(as.character(alldata$count_mismatch))),1,sum,na.rm=TRUE)
  othercols <- matrix(unlist(sapply(as.character(alldata$names_match),function(x){strsplit(x,split="_")})),nrow=nrow(alldata),byrow=TRUE)
  alldata <- cbind(othercols,total_count <- as.numeric(alldata$totalcount))
  mode(alldata) <-"numeric"
  alldata
}
