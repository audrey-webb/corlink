globalVariables(c("log_linear_match","log_linear_match_c","log_linear_mismatch","log_linear_mismatch_c"))
print_simplelmform <- function(L) {
        form <- c("~")
        for (lower in 1:L) {
            form <- paste(form, "var",lower, "+", sep = "")
        }
          for (lower in 1:(L-2)) {
            for(upper in (lower+1):L)
            form <- paste(form, "var",lower, "*","var",upper, "+", sep = "")
        }
        form <- paste(form, "var",(L-1),"*","var",L , sep = "")
        form
}

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


EM_match_modelsearch <- function(count_dframe,m_1,u_1,m_2,u_2,fixedcol=c(2:9),p_init=0.5,tol=10^-5,maxit=10000,allterms=allterms){
  N <- sum(count_dframe[ncol(count_dframe)])
  L <- ncol(count_dframe)-1
  colnames(count_dframe) <- c(paste("var",1:L,sep=""),"counts")

  stuff <- EM_match_independence_v2(count_dframe,m_1,u_1,m_2,u_2,p_init=0.5,tol=10^-5,fixedcol=fixedcol)         ##  get starting values.
  p_new <- stuff$p
  p_old <- p_new
  gamma_current <- stuff$probs
  counter=1
  probs_match_new <- count_dframe$counts*gamma_current/sum(count_dframe$counts*gamma_current)
  probs_match_old <- probs_match_new
  probs_mismatch_new <- count_dframe$counts*(1-gamma_current)/sum(count_dframe$counts*(1-gamma_current))
  probs_mismatch_old <- probs_mismatch_new

  while((counter<=maxit & (max(max(abs(probs_match_new-probs_match_old)),max(abs(probs_mismatch_new-probs_mismatch_old)),abs(p_new-p_old))  > tol))|counter<=2){
    #if(counter>1) browser()
    probs_match_old <-  probs_match_new
    probs_mismatch_old <-  probs_mismatch_new
    p_old <- p_new

    for(i in 1:nrow(count_dframe)){

      gamma_current[i] <- p_old*probs_match_old[i]/(p_old*probs_match_old[i]  + (1-p_old)*probs_mismatch_old[i])


    }
    gamma_current[is.na(gamma_current)] <- 0
    p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)
    Ecounts_match = data.frame(cbind(count_dframe[,1:L],"counts"=count_dframe$counts*gamma_current))
    Ecounts_mismatch = data.frame(cbind(count_dframe[,1:L],"counts"=count_dframe$counts*(1-gamma_current)))
    #      Ecounts_mismatch$counts <- round(Ecounts_mismatch$counts,0)
    #  Ecounts_match$counts <- round(Ecounts_match$counts,0)

    ##  only do the following on first iteration:
    if(counter==1){
      termstoadd_match <- ""
      newterm=""
      logN <- log(sum(count_dframe$counts))
      fullterms =  print_simplelmform(L)
      rm(log_linear_match,pos=1)
      eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
      termsnotadded = allterms
      termsadded = c()
      while(newterm!="<none>"){
        eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
        #stuff <- add1(log_linear_match, uppermatch$formula)    #strange behvious
        BICvec <- numeric(length(termsnotadded)+1)
        BICvec <- numeric(length(termsnotadded)+1)
        BICvec[1] <- log_linear_match$deviance
        for(j in 1:length(termsnotadded)){
          eval(parse(text=paste("log_linear_match_c <- glm(counts ~ .",termstoadd_match,"+",termsnotadded[j],",poisson , Ecounts_match)",sep="")))
          BICvec[1+j] <- log_linear_match_c$deviance + logN*(log_linear_match$df.residual - log_linear_match_c$df.residual)
        }
        newterm = termsnotadded[which.min(BICvec[2:length(BICvec)])]
        termstoadd_match <- paste(termstoadd_match,"+",newterm,sep="")
        termsnotadded = setdiff(termsnotadded,newterm)
        termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,L=ncol(count_dframe)-1))
        if(BICvec[1] == min(BICvec)) newterm="<none>"

      }
       termstoadd_mismatch <- ""
      newterm=""
      logN <- log(sum(count_dframe$counts))
      fullterms =  print_simplelmform(L)
      rm(log_linear_mismatch,pos=1)
      eval(parse(text=paste("log_linear_mismatch <- glm(counts ~ .",termstoadd_mismatch,",poisson , Ecounts_mismatch)",sep="")))
      termsnotadded = allterms
      termsadded = c()
      while(newterm!="<none>"){
        eval(parse(text=paste("log_linear_mismatch <- glm(counts ~ .",termstoadd_mismatch,",poisson , Ecounts_mismatch)",sep="")))
        #stuff <- add1(log_linear_mismatch, uppermismatch$formula)    #strange behvious
        BICvec <- numeric(length(termsnotadded)+1)
        BICvec[1] <- log_linear_mismatch$deviance
        for(j in 1:length(termsnotadded)){
          eval(parse(text=paste("log_linear_mismatch_c <- glm(counts ~ .",termstoadd_mismatch,"+",termsnotadded[j],",poisson , Ecounts_mismatch)",sep="")))
          BICvec[1+j] <- log_linear_mismatch_c$deviance + logN*(log_linear_mismatch$df.residual-log_linear_mismatch_c$df.residual)
        }
        newterm = termsnotadded[which.min(BICvec[2:length(BICvec)])]
        termstoadd_mismatch <- paste(termstoadd_mismatch,"+",newterm,sep="")
        termsnotadded = setdiff(termsnotadded,newterm)
        termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,L=ncol(count_dframe)-1))            ###  add possible 3 way interactions
        if(BICvec[1] == min(BICvec)) newterm="<none>"

      }
    }
    if(counter > 1){
      eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
      eval(parse(text=paste("log_linear_mismatch <- glm(counts ~ .",termstoadd_mismatch,",poisson , Ecounts_mismatch)",sep="")))
    }

    probs_match_new <- numeric(2^L)
    probs_mismatch_new <- numeric(2^L)
    match_probs <-  log_linear_match$fitted/sum(log_linear_match$fitted)
    probs_match_new[as.numeric(names(match_probs))] <- match_probs
    mismatch_probs <-  log_linear_mismatch$fitted/sum(log_linear_mismatch$fitted)
    probs_mismatch_new[as.numeric(names(mismatch_probs))] <- mismatch_probs
    # probs_match_new[ probs_match_new < 10^-10] <- 10^-10
    # probs_mismatch_new[probs_mismatch_new < 10^-10] <- 10^-10
    counter=counter + 1
  }
  return(list(p=p_new, probs = gamma_current, model_match = log_linear_match, model_mismatch=log_linear_mismatch))
}


EM_match_independence_v2 <- function(count_dframe,m_1,u_1,m_2,u_2,fixedcol=c(2:9),p_init=0.5,tol=10^-5){
  N <- sum(count_dframe[ncol(count_dframe)])
  L <- ncol(count_dframe)-1
  colnames(count_dframe) <- c(paste("var",1:L,sep=""),"counts")

  gamma_current <- numeric(nrow(count_dframe))
  count_pattern_match <- numeric(nrow(count_dframe))
  count_pattern_nomatch <- numeric(nrow(count_dframe))


  #  start off the algorithm
  p_old <- p_init
  m1_old <- m_1
  m2_old <- m_2
  u1_old <- u_1
  u2_old <- u_2

  for(i in 1:nrow(count_dframe)){

    cols_missing <- (1:L)[count_dframe[i,1:L]=="2"]
    n_missing <- length(cols_missing)
    cols_matching <- (1:L)[count_dframe[i,1:L]=="1"]
    n_match <- length(cols_matching)
    cols_notmatching <- (1:L)[count_dframe[i,1:L]=="0"]
    n_notmatch <- length(cols_notmatching)

    gamma_current[i] <- p_init*prod(m2_old[cols_missing])*prod(m1_old[cols_matching])*prod((1-m1_old-m2_old)[cols_notmatching])/(p_init*prod(m2_old[cols_missing])*prod(m1_old[cols_matching])*prod((1-m1_old-m2_old)[cols_notmatching])+(1-p_init)*prod(u2_old[cols_missing])*prod(u1_old[cols_matching])*prod((1-u1_old-u2_old)[cols_notmatching]))


  }


  p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)

  #create data frames for counts conditionally for both matches and non matches
  probs_match = pmax(count_dframe$counts*gamma_current,10^-10)/sum(pmax(count_dframe$counts*gamma_current,10^-10))
  probs_mismatch = pmax(count_dframe$counts*(1-gamma_current),10^-10)/sum(pmax(count_dframe$counts*(1-gamma_current),10^-10))

  #run log linear model
  m1_new <- numeric(L)
  u1_new <- u1_old
  m2_new <- numeric(L)
  u2_new <- numeric(L)
  for(j in 1:L){
    m1_new[j] <- sum(as.numeric(count_dframe[,j]==1)*probs_match)
    m2_new[j] <- sum(as.numeric(count_dframe[,j]==2)*probs_match)
    if(!(j %in% fixedcol)) u1_new[j] <- sum(as.numeric(count_dframe[,j]==1)*probs_mismatch)     #####  Update first name probabilities as these consitute the blocks and would be poorly estimated otherwise
    if(!(j %in% fixedcol)) u2_new[j] <- sum(as.numeric(count_dframe[,j]==2)*probs_mismatch)
  }
  while(max(max(abs(u1_new-u1_old)),max(abs(m1_new-m1_old)),abs(p_new-p_old),max(abs(m2_new-m2_old)),max(abs(u2_new-u2_old)))  > tol){
    m1_old <-  m1_new
    u1_old <-  u1_new
    m2_old <-  m2_new
    u2_old <-  u2_new
    p_old <- p_new


    for(i in 1:nrow(count_dframe)){

      cols_missing <- (1:L)[count_dframe[i,1:L]=="2"]
      n_missing <- length(cols_missing)
      cols_matching <- (1:L)[count_dframe[i,1:L]=="1"]
      n_match <- length(cols_matching)
      cols_notmatching <- (1:L)[count_dframe[i,1:L]=="0"]
      n_notmatch <- length(cols_notmatching)


      gamma_current[i] <- p_old*prod(m2_old[cols_missing])*prod(m1_old[cols_matching])*prod((1-m1_old-m2_old)[cols_notmatching])/(p_old*prod(m2_old[cols_missing])*prod(m1_old[cols_matching])*prod((1-m1_old-m2_old)[cols_notmatching])+(1-p_old)*prod(u2_old[cols_missing])*prod(u1_old[cols_matching])*prod((1-u1_old-u2_old)[cols_notmatching]))


    }
    p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)
    for(j in 1:L){
      m1_new[j] <- sum(as.numeric(count_dframe[,j]==1)*probs_match)
      m2_new[j] <- sum(as.numeric(count_dframe[,j]==2)*probs_match)
      if(!(j %in% fixedcol)) u1_new[j] <- sum(as.numeric(count_dframe[,j]==1)*probs_mismatch)     #####  Update first name probabilities as these consitute the blocks and would be poorly estimated otherwise
      if(!(j %in% fixedcol)) u2_new[j] <- sum(as.numeric(count_dframe[,j]==2)*probs_mismatch)
    }
    probs_match = pmax(count_dframe$counts*gamma_current,10^-10)/sum(pmax(count_dframe$counts*gamma_current,10^-10))
    probs_mismatch = pmax(count_dframe$counts*(1-gamma_current),10^-10)/sum(pmax(count_dframe$counts*(1-gamma_current),10^-10))

  }
  return(list(p=p_new, probs = gamma_current,m1=m1_new,u1=u1_new))
}


EM_match_modelsearch_iu <- function(count_dframe,m_1,u_1,m_2,u_2,fixedcol=c(2:9),p_init=0.5,tol=10^-5,maxit=10000,allterms=allterms){
    N <- sum(count_dframe[ncol(count_dframe)])
    L <- ncol(count_dframe)-1
    colnames(count_dframe) <- c(paste("var",1:L,sep=""),"counts")

     stuff <- EM_match_independence_v2(count_dframe,m_1,u_1,m_2,u_2,p_init=0.5,tol=10^-5,fixedcol=fixedcol)         ##  get starting values.
     p_new <- stuff$p
     p_old <- p_new
     gamma_current <- stuff$probs
     counter=1
     probs_match_new <- count_dframe$counts*gamma_current/sum(count_dframe$counts*gamma_current)
     probs_match_old <- probs_match_new
      probs_mismatch_new <- count_dframe$counts*(1-gamma_current)/sum(count_dframe$counts*(1-gamma_current))
     probs_mismatch_old <- probs_mismatch_new

        while((counter<=maxit & (max(max(abs(probs_match_new-probs_match_old)),max(abs(probs_mismatch_new-probs_mismatch_old)),abs(p_new-p_old))  > tol))|counter<=2){
             flush.console()
                       #if(counter>1) browser()
                probs_match_old <-  probs_match_new
               probs_mismatch_old <-  probs_mismatch_new
              p_old <- p_new

              for(i in 1:nrow(count_dframe)){

              gamma_current[i] <- p_old*probs_match_old[i]/(p_old*probs_match_old[i]  + (1-p_old)*probs_mismatch_old[i])


          }
          gamma_current[is.na(gamma_current)] <- 0
        p_new <- sum(gamma_current*count_dframe$counts)/sum(count_dframe$counts)
           Ecounts_match = data.frame(cbind(count_dframe[,1:L],"counts"=count_dframe$counts*gamma_current))
            Ecounts_mismatch = data.frame(cbind(count_dframe[,1:L],"counts"=count_dframe$counts*(1-gamma_current)))
           #      Ecounts_mismatch$counts <- round(Ecounts_mismatch$counts,0)
          #  Ecounts_match$counts <- round(Ecounts_match$counts,0)

          ##  only do the following on first iteration:
       if(counter==1){
           termstoadd_match <- ""
            newterm=""
            logN <- log(sum(count_dframe$counts))
            fullterms =  print_simplelmform(L)
             rm(log_linear_match,pos=1)
            eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
          termsnotadded = allterms
            termsadded = c()
            while(newterm!="<none>"){
                         eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
                     #stuff <- add1(log_linear_match, uppermatch$formula)    #strange behvious
              BICvec <- numeric(length(termsnotadded)+1)
                BICvec <- numeric(length(termsnotadded)+1)
                   BICvec[1] <- log_linear_match$deviance
              for(j in 1:length(termsnotadded)){
                     eval(parse(text=paste("log_linear_match_c <- glm(counts ~ .",termstoadd_match,"+",termsnotadded[j],",poisson , Ecounts_match)",sep="")))
                     BICvec[1+j] <- log_linear_match_c$deviance + logN*(log_linear_match$df.residual - log_linear_match_c$df.residual)
              }
              newterm = termsnotadded[which.min(BICvec[2:length(BICvec)])]
                termstoadd_match <- paste(termstoadd_match,"+",newterm,sep="")
                termsnotadded = setdiff(termsnotadded,newterm)
                 termsnotadded <- c(termsnotadded,getnewindexes_all(newterm,L=ncol(count_dframe)-1))
               if(BICvec[1] == min(BICvec)) newterm="<none>"
                      flush.console()
                 }
             termstoadd_mismatch <- ""
            eval(parse(text=paste("log_linear_mismatch <- glm(counts ~ .",termstoadd_mismatch,",poisson , Ecounts_mismatch)",sep="")))
            termstoadd_mismatch <- "var1"
            if(L > 1){newterm=paste("var",2:L,sep="")
              for(i in 1:length(newterm)) termstoadd_mismatch <- paste(termstoadd_mismatch,"+",newterm[i],sep="")
           }
        }
        if(counter > 1){
              eval(parse(text=paste("log_linear_match <- glm(counts ~ .",termstoadd_match,",poisson , Ecounts_match)",sep="")))
             eval(parse(text=paste("log_linear_mismatch <- glm(counts ~ ",termstoadd_mismatch,",poisson , Ecounts_mismatch)",sep="")))
        }

        probs_match_new <- numeric(2^L)
        probs_mismatch_new <- numeric(2^L)
        match_probs <-  log_linear_match$fitted/sum(log_linear_match$fitted)
         probs_match_new[as.numeric(names(match_probs))] <- match_probs
        mismatch_probs <-  log_linear_mismatch$fitted/sum(log_linear_mismatch$fitted)
          probs_mismatch_new[as.numeric(names(mismatch_probs))] <- mismatch_probs
        # probs_match_new[ probs_match_new < 10^-10] <- 10^-10
       # probs_mismatch_new[probs_mismatch_new < 10^-10] <- 10^-10
        counter=counter + 1
        }
      return(list(p=p_new, probs = gamma_current, model_match = log_linear_match, model_mismatch=log_linear_mismatch))
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
