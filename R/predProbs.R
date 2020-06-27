# An R function for computing predicted probabilities for binary,
# ordinal, multinomial, & generalized ordinal logit & probit models
# Version 1.0.1, Tim Liao, University of Illinois, March 2020

predProbs <- function(model,specs,method="logit") {
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
  out1 <- matrix(NA,nrow=nrow(specs),ncol=n.int)
  out2 <- matrix(NA,nrow=nrow(specs),ncol=n.int+1)
  if (method=="mlogit") lps <- rep(NA,n.int)
  for (i in 1:nrow(specs)) {
    for (j in 1:n.int) {
      if ({method=="logit" | method=="probit"} & n.int>1)
        lp <- unlist(int[j]-sum(pam*specs[i,]))
      if ({method=="logit" | method=="probit"} & n.int==1)
        lp <- unlist(int[j]+sum(pam*specs[i,]))
      if (method=="logit")
        out1[i,j] <- exp(lp)/{1+exp(lp)}
      if (method=="probit")
        out1[i,j] <- pnorm(lp)
      if (method=="gologit") {
        lp <- unlist(int[j]+sum(pam[j,]*specs[i,]))
        out1[i,j] <- exp(lp)/{1+exp(lp)}
      }
      if (method=="mlogit") {
        for (k in 1:n.int)
        lps[k] <- unlist(int[k]+sum(pam[k,]*specs[i,]))
        lps <- sum(exp(lps))
        lp <- unlist(int[j]+sum(pam[j,]*specs[i,]))
        out1[i,j] <- exp(lp)/{1+lps}
      }
    }
  }
  out1 <- cbind(out1,rep(1,nrow(specs)))
  for (i in 1:{n.int+1}) {
    if (i==1) out2[,i] <- out1[,i]
    if (i>1 & method!="mlogit") out2[,i] <- out1[,i] - out1[,i-1]
    if (i>1 & i<{n.int+1} & method=="mlogit") out2[,i] <- out1[,i]
    if (i=={n.int+1} & method=="mlogit")
      out2[,i] <- out1[,i] - apply(out1[,1:{i-1}],1,sum)
  }
  cn <- rep("p",n.int+1)
  sn <- 1:{n.int+1}
  colName <- paste(cn,sn,sep="")
  if (n.int==1) colName <- c("p1","p0")
  if (method=="mlogit")
    colnames(out2) <- c(colName[1:{ncol(out2)-1}],"pRef")
  else colnames(out2) <- colName[1:ncol(out2)]
  return(data.frame(out2))
}
