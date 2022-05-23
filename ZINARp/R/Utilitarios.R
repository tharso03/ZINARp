
#' @importFrom stats dbinom
#' @importFrom stats dpois
#' @importFrom stats rbeta
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @importFrom stats rgamma
#' @importFrom stats rpois

# Innovation
####################

genINARp <- function(alpha,lambda,n)
{
  p <- length(alpha)
  y0 <- rep(lambda/(1-sum(alpha)),p)
  burn <- 5*n
  y <- numeric()
  y[1:p] <- round(y0,0)
  func_opbin <- function(contador)
  {
    x_prob = lista$prob[contador]
    x_y = lista$y[contador]
    rbinom(1,size=x_y,prob=x_prob)
  }
  k = p +1
  while(k <= (n+burn))
  {
    lista = list(y=y[(k-1):(k-p)],prob = alpha)
    opbin = sum(sapply(1:p,func_opbin))
    zip = rpois(1,lambda)
    y[k] = opbin + zip
    k = k + 1
  }
  y[(burn + 1):(n+burn)]
}

simuZIP <- function(n,pii,lambda)
{
  y <- c()
  U <- runif(n)
  for(i in 1:n)
  {
    y[i]<- ifelse(pii>=U[i],0,rpois(1,lambda))
  }
  return(y)
}

genZINARp <- function(alpha,pii,lambda,n)
{
  p <- length(alpha)
  y0 <- rep(lambda/(1-sum(alpha)),p)
  burn <- 5*n
  y <- numeric()
  y[1:p] <- round(y0,0)
  k = p +1
  while(k <= (n+burn))
  {
    opbin <- c(0)
    for(i in 1:p)
    {
      opbin[i] <- rbinom(1,y[k-i], alpha[i])
    }
    zip  = simuZIP(1,pii,lambda)
    y[k] = sum(opbin) + zip
    k = k + 1
  }
  y1 <- y[(burn + 1):(n+burn)]
  return(y1)
}


# lambda / y,alpha - Inarp
##########################

MH.StI<- function (y,lastSt,lastZt,alpha,lambda,t,p=p)
{
  n <- length(y)
  StCand <- matrix(0,1,p)
  aux <- 0
  while(aux==0)
  {
    for(k in 1:p)
    {
      if(y[t-k]==0)
      {
        StCand[1,k] = 0
      }
      else
      {
        StCand[1,k] <- rbinom(1,y[t-k],alpha[k])
      }
    }
    if(y[t]<sum(StCand))
    {
      aux<- 0
    }
    else
    {
      ZtCand <- y[t]-sum(StCand)
      aux <- 1
    }
  }
  unif <- runif(1)
  rej <- 0
  if(unif*(lambda^(lastZt)/factorial(lastZt))< (lambda^(ZtCand)/factorial(ZtCand)))
  {
    next.St <- StCand
    next.Zt <- ZtCand
  }
  else{
    next.St <- lastSt
    next.Zt <- lastZt
    rej <- 1
  }
  return(list(St1=next.St,Zt1=next.Zt,rejeitou=rej))
}


# Densidade do ZIP(rho,lambda)
##############################

DensZIP <- function(y,rho,lambda)
{
  if(y==0){
    dens <- rho + (1-rho)*exp(-lambda)}
  else{
    dens <- (1-rho)*dpois(y,lambda=lambda)}
  return(dens)
}


# St / x,alpha, lambda,rho,wt
##########################

MH.St <- function (y,lastSt,lastZt,lastWt,alpha,lambda,rho,t,p=p)
{
  n <- length(y)
  StCand <- matrix(0,1,p)
  aux <- 0
  while(aux==0)
  {
    for(k in 1:p)
    {
      if(y[t-k]==0)
      {
        StCand[1,k] = 0
      }
      else
      {
        StCand[1,k] <- rbinom(1,y[t-k],alpha[k])
      }
    }
    if(y[t]<sum(StCand))
    {
      aux<- 0
    }
    else
    {
      ZtCand <- y[t]-sum(StCand)
      WtCand <- GerWt(ZtCand,lambda=lambda,rho=rho)
      aux <- 1
    }
  }
  unif <- runif(1)
  rej <- 0
  if(unif*((lambda)^(lastZt*(1-lastWt)))*dbinom(lastWt,1,rho)*DensWt(WtCand,ZtCand,lambda,rho)*exp(-lambda*(1-lastWt))/factorial((lastZt))< ((lambda)^(ZtCand*(1-WtCand)))*dbinom(WtCand,1,rho)*DensWt(lastWt,lastZt,lambda,rho)*exp(-lambda*(1-WtCand))/factorial((ZtCand)))
  {
    next.St <- StCand
    next.Zt <- ZtCand
    next.Wt <- WtCand
  }
  else{
    next.St <- lastSt
    next.Zt <- lastZt
    next.Wt <- lastWt
    rej <- 1
  }
  return(list(St1=next.St,Zt1=next.Zt,Wt1=next.Wt,rejeitou=rej))
}

# Wt / Zt, lambda,rho,x
##########################

GerWt <- function(Zt,lambda,rho)
{
  n1 <- length(Zt)
  wt <- c()
  p1 <- rho/(rho + (1-rho)*exp(-lambda))
  for(i in 1:n1)
  {
    if(Zt[i]>0)
    {
      wt[i] <- 0
    }
    else
    {
      wt[i] <- rbinom(1,1,p1)
    }
  }
  return(wt)
}

DensWt <- function(Wt,Zt,lambda,rho)
{
  n1 <- length(Zt)
  dens <- c()
  p1 <- rho/(rho + (1-rho)*exp(-lambda))
  for(i in 1:n1)
  {
    if(Zt[i]>0)
    {
      dens[i] <- ifelse(Wt[i]==0,1,0)
    }
    else
    {
      dens[i] <- ifelse(Wt[i]==1,dbinom(1,1,p1),dbinom(0,1,p1))
    }
  }
  return(dens)
}

DensWtZt <- function(Wt,Zt,lambda)     ## Densidade de Zt / Wt
{
  n1 <- length(Zt)
  dens <- c()
  if(Wt==1)
  {
    dens <- ifelse(Zt==0,1,0)
  }
  else
  {
    dens <- dpois(Zt,lambda)
  }
  return(dens)
}

GerZtWt <- function(Wt,lambda)
{
  Zt <- c()
  if(Wt==1)
  {
    Zt <- 0
  }
  else
  {
    Zt <- rpois(1,lambda)
  }
  return(Zt)
}

## Gibbs sampler

GibbsINARp <-function(x,p,iter=1000,thin=10,burn=0.2)
{
  x <- as.vector(x)
  n <- length(x)
  pb <- progress::progress_bar$new(total=iter)

  ###### Hiper par?metros: comum ######
  A <- 0.001
  B <- 0.001

  ###### Matriz de resultados ######
  alphaG <- matrix(0,iter,p)
  lambdaG <- matrix(0,iter,1)
  St <- matrix(0,n,p)
  Zt <- matrix(0,n,1)

  ###### Valores Iniciais ######

  alphaG[1,] <- rep(0.1,p)
  lambdaG[1] <- runif(1,0,5)
  St[(p+1),] <- rep(0,p)
  Zt[p+1] <- x[p+1]- sum(St[(p+1),])
  for (j in 2:iter)
  {
    i <- p+1
    while(i <=n)
    {
      Aux <- MH.StI(x,lastSt=St[i,],lastZt=Zt[i],alpha=alphaG[j-1,],lambda=lambdaG[j-1],t=i,p=p)
      St[i,] <- Aux$St1
      Zt[i] <- Aux$Zt1
      i <- i+1
    }
    Aux2 <- matrix(0,n,p)
    aux1 = 0
    while(aux1==0)
    {
      for(l in 1:p)
      {
        for(k in (p+1):n){Aux2[k,l] <-x[k-l]-St[k,l]}
        alphaG[j,l] <- rbeta(1,sum(St[(p+1):n,l])+1,sum(Aux2[,l])+1)
      }
      aux1 <- ifelse(sum(alphaG[j,])<1,1,0)
    }
    lambdaG[j] <- rgamma(1,shape=sum(Zt[(p+1):n])+A,rate=(n-p+B))
    pb$tick()
    S1 <- c()
  }
  ###burning: 20% - Default
  initial<-ceiling(iter*burn)
  set<-seq(initial,iter,thin)
  return(list(alpha=alphaG[set,],lambda=lambdaG[set]))
}


GibbsZINARp <-function(x,p,iter=1000,thin=10,burn=0.2)
{
  x <- as.vector(x)
  n <- length(x)
  pb <- progress::progress_bar$new(total=iter)

  ###### Hiper par?metros: comum ######
  A <- 0.001
  B <- 0.001

  ###### Matriz de resultados ######
  alphaG <- matrix(0,iter,p)
  lambdaG <- matrix(0,iter,1)
  rhoG <- matrix(0,iter,1)
  St <- matrix(0,n,p)
  Zt <- matrix(0,n,1)
  Wt <- matrix(0,n,1)



  ###### Valores Iniciais ######

  alphaG[1,] <- rep(0.1,p)
  lambdaG[1] <- runif(1,0,5)
  rhoG[1] <- runif(1)
  St[(p+1),] <- rep(0,p)
  Zt[p+1] <- x[p+1]- sum(St[(p+1),])
  Wt[p+1] <- ifelse(Zt[p+1]==0,rbinom(1,1,rhoG[1]),0)
  for (j in 2:iter)
  {
    i <- p+1
    while(i <=n)
    {
      Aux <- MH.St(y=x,lastSt=St[i,],lastZt=Zt[i],lastWt=Wt[i],alpha=alphaG[j-1,],lambda=lambdaG[j-1],rho=rhoG[j-1],t=i,p=p)
      St[i,] <- Aux$St1
      Zt[i] <- Aux$Zt1
      Wt[i] <- Aux$Wt1
      i <- i+1
    }
    Aux2 <- matrix(0,n,p)
    aux1 = 0
    while(aux1==0)
    {
      for(l in 1:p)
      {
        for(k in (p+1):n){Aux2[k,l] <-x[k-l]-St[k,l]}
        alphaG[j,l] <- rbeta(1,sum(St[(p+1):n,l])+1,sum(Aux2[,l])+1)
      }
      aux1 <- ifelse(sum(alphaG[j,])<1,1,0)
    }
    lambdaG[j] <- rgamma(1,shape=sum(Zt[(p+1):n]*(1-Wt[(p+1):n]))+A,rate=sum(1-Wt[(p+1):n])+B)
    rhoG[j] <- rbeta(1,shape1=sum(Wt[(p+1):n])+1,shape2=sum(1-Wt[(p+1):n])+1)
    S1 <- c()
    pb$tick()
  }
  ###burning: 20%
  initial<-ceiling(iter*burn)
  set<-seq(initial,iter,thin)
  return(list(alpha=alphaG[set,],lambda=lambdaG[set],rho=rhoG[set]))
}



