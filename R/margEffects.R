# An R function for computing predicted probabilities for binary,
# ordinal, multinomial, & generalized ordinal logit & probit models
# Version 1.0.1, Tim Liao, University of Illinois, March 2020

margEffects <- function(model,specs,effect=1,method="logit") {
  if (method=="logit" | method=="probit") {
    n.int <- length(unique(model$y)) - 1
    n.ind <- length(model$coef) - n.int
    int <- model$coef[1:n.int]
    pam <- model$coef[{n.int+1}:{n.int+n.ind}]
  }
  if (method=="gologit" | method=="mlogit") {
    n.int <- ncol(predict(model))
    cf <- coef(model)
    n.ind <- {length(cf)-n.int}/n.int
    int <- coef(model)[1:n.int]
    pam <- matrix(NA,nrow=n.int,ncol=n.ind)
    t <- rep(NA,n.ind)
    for (i in 1:n.ind) t[i] <- rep(0+i*n.int)
    for (i in 1:n.int)
      pam[i,] <- cf[i+t]
  }
  if (method=="logit" | method=="probit" | method=="gologit")
    den=prob <- matrix(NA,nrow=nrow(specs),ncol=n.int)
  if (method=="mlogit")
    prob <- matrix(NA,nrow=nrow(specs),ncol=n.int+1)
  if (method=="mlogit") lps <- rep(NA,n.int)
  for (i in 1:nrow(specs)) {
    for (j in 1:n.int) {
      if ({method=="logit" | method=="probit"} & n.int>1)
        lp <- unlist(int[j]-sum(pam*specs[i,]))
      if ({method=="logit" | method=="probit"} & n.int==1)
        lp <- unlist(int[j]+sum(pam*specs[i,]))
      if (method=="gologit")
        lp <- unlist(int[j]+sum(pam[j,]*specs[i,]))
      if (method=="logit" | method=="gologit")
        prob[i,j] <- exp(lp)/{1+exp(lp)}
      if (method!="mlogit")
        den[i,j] <- prob[i,j]*(1-prob[i,j])
      if (method=="probit")
        den[i,j] <- dnorm(lp)
      if (method=="mlogit") {
        for (k in 1:n.int)
          lps[k] <- unlist(int[k]+sum(pam[k,]*specs[i,]))
        lps <- sum(exp(lps))
        lp <- unlist(int[j]+sum(pam[j,]*specs[i,]))
        if (j<=n.int)
          prob[i,j] <- exp(lp)/{1+lps}
        if (j==n.int)
          prob[i,j+1] <- 1/{1+lps}
      }
    }
  }
  out <- matrix(NA,nrow=nrow(specs),ncol=n.int+1)
  if (method=="logit" | method=="probit") {
    for (i in 1:{n.int+1}) {
      if (n.int>1) {
        if (i==1) out[,i] <- -den[,i]*pam[effect]
        else if (i>1 & i<{n.int+1})
          out[,i] <- {den[,i-1] - den[,i]}*pam[effect]
        else if (i=={n.int+1})
          out[,i] <- den[,i-1]*pam[effect]
      }
      if (n.int==1) {
        if (i==1) out[,i] <- den[,i]*pam[effect]
        if (i>1) out[,i] <- -den[,i-1]*pam[effect]
      }
    }
  }
  if (method=="gologit") {
    for (i in 1:{n.int+1})
      if (i==1) out[,i] <- den[,i]*pam[i,effect]
      else if (i>1 & i<{n.int+1}) {
        out[,i] <- -den[,i-1]*pam[i-1,effect] + den[,i]*pam[i,effect] }
      else if (i=={n.int+1})
        out[,i] <- -den[,i-1]*pam[i-1,effect]
  }
  si <- matrix(NA,nrow=nrow(specs),ncol=n.int)
    if (method=="mlogit") {
    for (i in 1:n.int)
      si[,i] <- prob[,i]*pam[i,effect]
    s <- apply(si,1,sum)
    for (j in 1:n.int) {
        out[,j] <- prob[,j]*(pam[j,effect]-s)
        if (j==n.int)
          out[,j+1] <- prob[,j+1]*(0-s)
    }
  }

  cn <- rep("p",n.int+1)
  sn <- 1:{n.int+1}
  colName <- paste(cn,sn,sep="")
  if (n.int==1) colName <- c("p1","p0")
  if (method=="mlogit")
    colnames(out) <- c(colName[1:{ncol(out)-1}],"pRef")
  if (method!="mlogit") colnames(out) <- colName[1:ncol(out)]

  return(data.frame(out))
}
