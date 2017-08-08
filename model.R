kde <- function(x, w = NA){
  if(sum(is.na(w)) >0){
    w <- seq(1,1,length.out = length(x))
  }
  n <- 100
  r <- hist(x,n,plot = F)
  p <- ceiling((x - min(r$breaks)) / (r$breaks[2] - r$breaks[1]))
  mids <- r$mids
  p[p<1] <- 1
  p[p>length(mids)] <-length(mids)
  weights <- seq(0,0,length.out = length(mids))
  for(i in 1:length(mids)){
    weights[i] <- sum(w[p==i])
  }
  m <- sum(weights * mids) / sum(weights)
  s <- sqrt(sum((mids - m)^2) / sum(weights))
  h <- 0.9 * min(s, quantile(mids)[4] - quantile(mids)[2])* length(mids)^(-0.2)
  
  dist <- matrix(rep(mids,length(mids)),nrow=length(mids)) - matrix(rep(mids,length(mids)),nrow=length(mids),byrow = T)
  dens <- dnorm(dist / h)
  y <- apply(dens * matrix(rep(weights,length(mids)),nrow=length(mids),byrow=T), 1, sum) / h / sum(weights)
  
  y <- y[p]
  return(y)
}

sLDA <- function(z,v){
  n <- length(z)
  mu1 <- apply(v * z, 2, sum) / sum(z)
  mu0 <- apply(v * (1-z), 2, sum) / sum(1-z)
  mu <- (sum(z) * mu1 + sum(1 - z) * mu0) / n
  Sb <- 1/n*(matrix((mu0-mu),ncol=1) %*% matrix((mu0-mu),nrow=1)*sum(1-z) + matrix((mu1-mu),ncol=1) %*% matrix((mu1-mu),nrow=1) * sum(z))
  S <- 1/n*( t( z * (v - matrix(rep(mu1,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu1,n),nrow=n,byrow = T))  + t( (1-z) * (v - matrix(rep(mu0,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu0,n),nrow=n,byrow = T)) )
  r <- eigen(solve(S + 1e-9 * diag(dim(S)[1])) %*% Sb)
  lda.v <- r$vectors[,1]
  y <- as.numeric(v %*% lda.v)
  return(y)
}

sLDA.2 <- function(z,v){
  n <- length(z)
  mu1 <- apply(v * z, 2, sum) / sum(z)
  mu0 <- apply(v * (1-z), 2, sum) / sum(1-z)
  mu <- (sum(z) * mu1 + sum(1 - z) * mu0) / n
  Sb <- 1/n*(matrix((mu0-mu),ncol=1) %*% matrix((mu0-mu),nrow=1)*sum(1-z) + matrix((mu1-mu),ncol=1) %*% matrix((mu1-mu),nrow=1) * sum(z))
  S <- 1/n*( t( z * (v - matrix(rep(mu1,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu1,n),nrow=n,byrow = T))  + t( (1-z) * (v - matrix(rep(mu0,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu0,n),nrow=n,byrow = T)) )
  r <- eigen(solve(S + 1e-9 * diag(dim(S)[1])) %*% Sb)
  lda.v <- as.numeric(r$vectors[,1])
  return(lda.v)
}

EMinfer <- function(pvals){
  gene.num <- length(pvals)
  
  alpha0 <- 1
  alpha1 <- 0.2
  pi0 <- 0.6
  pi1 <- 0.4
  for(i in 1:200){
    tmp1 <- pi1*alpha1*pvals^(alpha1-1)
    tmp0 <- pi0*alpha0*(pvals)^(alpha0-1)
    proZ <- tmp1/(tmp1+tmp0)
    
    pi0 <- sum(1-proZ)/length(proZ)
    pi1 <- sum(proZ) / length(proZ)
    
    alpha0 <- - sum(1 - proZ)/ sum( (1-proZ)* log(pvals) )
    alpha1 <- - sum(proZ)/ sum( proZ* log(pvals) )
    
  }
  result <- list(alpha0 = alpha0, alpha1 = alpha1, pi = pi1, post=proZ)
  return(result)
}


model <- function(pvalues,v,verbose=F){
  n <- dim(v)[1]
  d <- dim(v)[2]
  if(length(pvalues) != n){
    print('Error : length of pvalues not equal to dim(v)[1] !')
    return
  }
  
  r <- EMinfer(pvalues)
  alpha0 <- r$alpha0
  alpha1 <- r$alpha1
  pi <- r$pi
  z <- r$post
  y <- sLDA(z,v)
  old.z <- z
  
  density0 <- seq(0,0,length.out = n)
  density1 <- seq(0,0,length.out = n)
  
  loglik.history <- c()
  i <- 0
  while(1){
    i <- i + 1
    
    density0 <- kde(y,1-z)
    density1 <- kde(y,z)
    
    # E-step
    tmp0 <- log(1-pi) + log(alpha0) + (alpha0-1)*log(pvalues) + log(density0)
    tmp1 <- log(pi) + log(alpha1) + (alpha1-1)*log(pvalues) + log(density1)
    z <- 1 / (1 + exp(tmp0 - tmp1))
    
    # M-step
    pi <- sum(z) / n
    alpha0 <- - sum(1-z) / sum((1-z) * log(pvalues))
    alpha1 <- - sum(z) / sum(z * log(pvalues))
    
    if(sum(abs(old.z - z)) < 1e-6) break;
    if(i > 1000) break;
    old.z <- z
  }
  
  result <- list(post=z, y = y, pi=pi, alpha0=alpha0, alpha1=alpha1, density0 = density0[order(y)], density1 = density1[order(y)])
  return(result)
}

model.M <- function(pvalues,v,verbose=F){
  n <- dim(v)[1]
  d <- dim(v)[2]
  if(length(pvalues) != n){
    print('Error : length of pvalues not equal to dim(v)[1] !')
    return
  }
  
  r <- EMinfer(pvalues)
  alpha0 <- r$alpha0
  alpha1 <- r$alpha1
  pi <- r$pi
  z <- r$post
  old.z <- z
  
  density0 <- matrix(0,nrow = n, ncol = d)
  density1 <- matrix(0,nrow = n, ncol = d)
  
  loglik.history <- c()
  i <- 0
  while(1){
    i <- i + 1
    
    for(k in 1:d){
      density0[,k] <- kde(v[,k],1-z)
      density1[,k] <- kde(v[,k],z)
    }
    
    # E-step
    tmp0 <- log(1-pi) + log(alpha0) + (alpha0-1)*log(pvalues) + apply(log(density0),1,sum)
    tmp1 <- log(pi) + log(alpha1) + (alpha1-1)*log(pvalues) + apply(log(density1),1,sum)
    z <- 1 / (1 + exp(tmp0 - tmp1))
    
    # M-step
    pi <- sum(z) / n
    alpha0 <- - sum(1-z) / sum((1-z) * log(pvalues))
    alpha1 <- - sum(z) / sum(z * log(pvalues))
    
    if(sum(abs(old.z - z)) < 1e-6) break;
    if(i > 1000) break;
    old.z <- z
  }
  
  result <- list(post=z, pi=pi, alpha0=alpha0, alpha1=alpha1)
  return(result)
}




model.LR <- function(pvalues,v,verbose=F, lambda=0, cv = 0){
  
  # initilization
  n <- dim(v)[1]
  d <- dim(v)[2]
  if(length(pvalues) != n){
    print('Error : length of pvalues not equal to dim(v)[1] !')
    return
  }
  
  r <- EMinfer(pvalues)
  w <- seq(0,0,length.out = d)
  alpha0 <- r$alpha0
  alpha1 <- r$alpha1
  b <- 0
  z <- r$post
  gamma <- 1
  
  pvalues[pvalues < 1e-50] <- 1e-50
  
  loglik.history <- c()
  i <- 0
  while(1){
    i <- i + 1
    # E-step
    prior <- 1/(1+exp(- v %*% w - b))
    prior <- prior[,1]
    
    tmp0 <- log(1 - prior) + log(alpha0) + (alpha0 - 1) * log(pvalues)
    tmp1 <- log(prior) + log(alpha1) + (alpha1 - 1) * log(pvalues)
    z <- 1 / (1 + exp(tmp0 - tmp1))
    
    max.tmp <- max(c(tmp1,tmp0))
    loglik <- sum(log(exp(tmp1-max.tmp) + exp(tmp0-max.tmp)) + max.tmp) - lambda/2 * n* sum(w*w)
    loglik.history <- c(loglik.history,loglik)
    
    # M-step
    alpha1 <- - sum(z) / sum(z * log(pvalues))
    alpha0 <- - sum(1-z) / sum((1-z) * log(pvalues))
    
    grad.w <- apply(v * (z - gamma * prior + (gamma -1) * z * prior),2,sum) - lambda * n * w 
    grad.b <- sum(z - gamma * prior + (gamma - 1) * z * prior)
    
    H <- matrix(0,ncol = d+1,nrow = d+1)
    tmp <- prior*(prior-1) #* idx
    H[1:d,1:d] <- t(v) %*% (v * (gamma * tmp + (1-gamma) * z * tmp) ) - lambda * n * diag(d)
    H[d+1,d+1] <- sum(gamma * tmp + (1-gamma) * z * tmp)
    H[1:d,d+1] <- apply((v * (gamma * tmp + (1-gamma)*z*tmp ) ),2,sum)
    H[d+1,1:d] <- apply((v * (gamma * tmp + (1-gamma)*z*tmp ) ),2,sum)
    
    H <- H + 1e-6 * diag(d+1)
    step <- solve(H, c(grad.w,grad.b))
    w <- w - step[1:d] 
    b <- b - step[d+1]
    
    if(verbose) print(sprintf('iter %02d , log(likelihood) %.6f',i,loglik))
    
    if(i > 10 && abs(loglik.history[i-10] - loglik.history[i]) < 1e-6 ) break;
    if(i >= 500) break;
    if(i > 1 && loglik.history[i] < loglik.history[i-1]) break;
    #print(roc(data[,'lab'],z)$auc[[1]])
  }
  
  result <- list(loglik=loglik.history, alpha0=alpha0,alpha1=alpha1,w=w,b=b,post=z)
  
  return(result)
}


model.NB <- function(pvalues,v,verbose=F){
  # initilization
  n <- dim(v)[1]
  d <- dim(v)[2]
  if(length(pvalues) != n){
    print('Error : length of pvalues not equal to dim(v)[1] !')
    return
  }
  
  
  r <- EMinfer(pvalues)
  alpha0 <- r$alpha0
  alpha1 <- r$alpha1
  pi <- r$pi
  mu0 <- seq(0,0,length.out = d)
  mu1 <- seq(0,0,length.out = d)
  sigma <- seq(1,1,length.out = d)
  z <- seq(0,0,length.out = n)
  
  pvalues[pvalues < 1e-50] <- 1e-50
  
  loglik.history <- c()
  i <- 0
  while(1){
    i <- i + 1
    
    # E-step
    tmp0 <- log(1-pi) + log(alpha0) + (alpha0-1)*log(pvalues) - 1/2 * (sum(log(sigma)) + apply((v - matrix(rep(mu0,n),nrow=n,byrow = T))^2/matrix(rep(sigma,n),nrow=n,byrow=T),1,sum) )
    tmp1 <- log(pi) + log(alpha1) + (alpha1-1)*log(pvalues) - 1/2 * (sum(log(sigma)) + apply((v - matrix(rep(mu1,n),nrow=n,byrow = T))^2/matrix(rep(sigma,n),nrow=n,byrow=T),1,sum) )
    z <- 1 / (1 + exp(tmp0 - tmp1))
    
    max.tmp <- max(c(tmp1,tmp0))
    loglik <- sum(log(exp(tmp1-max.tmp) + exp(tmp0-max.tmp)) + max.tmp)
    loglik.history <- c(loglik.history,loglik)
    
    # M-step
    pi <- sum(z) / n
    alpha0 <- - sum(1-z) / sum((1-z) * log(pvalues))
    alpha1 <- - sum(z) / sum(z * log(pvalues))
    
    mu0 <- apply(v * (1-z), 2, sum) / sum(1-z)
    mu1 <- apply(v * z, 2, sum) / sum(z)
    
    sigma <- apply(1/n*(z * (v - matrix(rep(mu1,n),nrow=n,byrow = T))^2 + (1-z) * (v - matrix(rep(mu0,n),nrow=n,byrow = T))^2), 2, sum)
    
    if(verbose) print(sprintf('iter %02d , log(likelihood) %.6f',i,loglik))
    
    if(i > 10 && abs(loglik.history[i-10] - loglik.history[i]) < 1e-6 ) break;
    if(i >= 5000) break;
    if(i > 1 && loglik.history[i] < loglik.history[i-1]) break;
    #print(roc(data[,'lab'],z)$auc)
  }
  
  result <- list(loglik=loglik.history, post=z, alpha0=alpha0,alpha1=alpha1,pi=pi,mu0=mu0,mu1=mu1,sigma=sigma)
  return(result)
}

model.MVN <- function(pvalues,v,verbose=F){
  # initilization
  n <- dim(v)[1]
  d <- dim(v)[2]
  if(length(pvalues) != n){
    print('Error : length of pvalues not equal to dim(v)[1] !')
    return
  }
  
  r <- EMinfer(pvalues)
  alpha0 <- r$alpha0
  alpha1 <- r$alpha1
  pi <- r$pi
  mu0 <- seq(0,0,length.out = d)
  mu1 <- seq(0,0,length.out = d)
  sigma <- diag(d)
  z <- seq(0,0,length.out = n)
  
  loglik.history <- c()
  i <- 0
  while(1){
    i <- i + 1
    
    # E-step
    tmp0 <- log(1-pi) + log(alpha0) + (alpha0-1)*log(pvalues) - 1/2 * (sum(log(det(sigma))) + apply( (v-matrix(rep(mu0,n),nrow=n,byrow = T)) %*% solve(sigma) * (v-matrix(rep(mu0,n),nrow=n,byrow = T)),1,sum) )
    tmp1 <- log(pi) + log(alpha1) + (alpha1-1)*log(pvalues) - 1/2 * (sum(log(det(sigma) )) + apply( (v-matrix(rep(mu1,n),nrow=n,byrow = T)) %*% solve(sigma) * (v-matrix(rep(mu1,n),nrow=n,byrow = T)),1,sum) )
    z <- 1 / (1 + exp(tmp0 - tmp1))
    
    
    max.tmp <- max(c(tmp1,tmp0))
    loglik <- sum(log(exp(tmp1-max.tmp) + exp(tmp0-max.tmp) + 1e-320) + max.tmp)
    loglik.history <- c(loglik.history,loglik)
    
    # M-step
    pi <- sum(z) / n
    alpha0 <- - sum(1-z) / sum((1-z) * log(pvalues))
    alpha1 <- - sum(z) / sum(z * log(pvalues))
    
    mu0 <- apply(v * (1-z), 2, sum) / sum(1-z)
    mu1 <- apply(v * z, 2, sum) / sum(z)
    
    sigma <- 1/n*( t( z * (v - matrix(rep(mu1,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu1,n),nrow=n,byrow = T))  + t( (1-z) * (v - matrix(rep(mu0,n),nrow=n,byrow = T))) %*% (v - matrix(rep(mu0,n),nrow=n,byrow = T)) )
    
    if(verbose) print(sprintf('iter %02d , log(likelihood) %.6f',i,loglik))
    
    if(i > 10 && abs(loglik.history[i-10] - loglik.history[i]) < 1e-6 ) break;
    if(i >= 5000) break;
    if(i > 1 && loglik.history[i] < loglik.history[i-1]) break;
    #print(roc(data[,'lab'],z)$auc)
  }
  
  result <- list(loglik=loglik.history, post=z, alpha0=alpha0,alpha1=alpha1,pi=pi,mu0=mu0,mu1=mu1,sigma=sigma)
  return(result)
}