# pkgdown::build_site()
# R CMD check --as-cran PRECAST_1.0.tar.gz
# devtools::check_win_release()
# iDR.SC <- function(...) UseMethod("iDR.SC")
model_set <- function(Sigma_equal=FALSE, Sigma_diag=TRUE,mix_prop_heter=TRUE,
                      error_heter=TRUE, Sp2=TRUE, wpca_int=FALSE,int.model='EEE',
                      coreNum = 1, coreNum_int=coreNum,
                      beta_grid=seq(0.2,4, by=0.2),
                      maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE, seed=1){
  para_settings <- list(Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                        mix_prop_heter=mix_prop_heter,
                        error_heter=error_heter, Sp2=Sp2, wpca_int=wpca_int,int.model=int.model,
                        coreNum = coreNum, coreNum_int=coreNum_int,
                        beta_grid= beta_grid,
                        maxIter_ICM=maxIter_ICM,maxIter= maxIter, epsLogLik=epsLogLik,
                        verbose=verbose, seed=seed)
  return(para_settings)
  
}
ICM.EM_structure  <- function(XList,  K, AdjList, q=15,parameterList=NULL){
  
  if(is.null(parameterList)){
    parameterList <- model_set()
  }
  ### initialize variables
  beta_grid <- maxIter_ICM<- maxIter<- epsLogLik<- verbose<-mix_prop_heter<- Sigma_equal<- Sigma_diag<- NULL
  error_heter<-Sp2<- wpca_int<-int.model<- seed<- coreNum <- coreNum_int <- NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }
  resList <- ICM.EM(XList, q, K, AdjList=AdjList, 
              beta_grid=beta_grid,
              maxIter_ICM=maxIter_ICM,maxIter=maxIter, epsLogLik=epsLogLik, verbose=verbose,
              mix_prop_heter=mix_prop_heter, Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
              error_heter=error_heter, Sp2=Sp2,
              wpca_int=wpca_int,int.model=int.model, seed=seed,coreNum = coreNum, coreNum_int=coreNum_int)
  return(resList)
}
ICM.EM <- function(XList, q, K, AdjList=NULL,  Adjlist_car=NULL, posList = NULL, platform = "ST",
                   beta_grid=seq(0.2,4, by=0.2),
                   maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE,
                   mix_prop_heter=TRUE, Sigma_equal=FALSE, Sigma_diag=TRUE,error_heter=TRUE, Sp2=TRUE,
                   wpca_int=FALSE,int.model='EEE', seed=1,coreNum = 1, coreNum_int=coreNum){
  
  
  
  
  
  
  r_max <- length(XList)
  
  pf <- unlist(lapply(XList, ncol))
  # if(is.null(posList) && is.null(AdjList)) stop("Check the values of arguments!")
  if(sum(pf!=pf[1]) > 0) stop("ncols of XList's components must be same!")
  if(is.null(posList) && is.null(AdjList)) stop("posList or AdjList must be provided one!")
  n <- sum(sapply(XList, nrow))
  p <-  pf[1]
  message("Intergrative data info.: ", r_max, ' samples, ', p , " genes X ", n, " spots------")
  message("PRECAST model setting: ", 'error_heter=',error_heter,", Sigma_equal=",Sigma_equal,
          ", Sigma_diag=", Sigma_diag, ', mix_prop_heter=', mix_prop_heter)
  
  if(is.null(AdjList) && !is.null(posList))
    AdjList <- lapply(posList, getAdj, platform=platform)
  if(any(sapply(AdjList, nrow)!=sapply(XList, nrow)))  
    stop('The dimension of Adj matrix does not match that of X!\n')
  
  if(is.null(Adjlist_car)) Adjlist_car <- AdjList
  #XList <- lapply(XList, scale, scale=scale_flag)
  Xmat <- NULL
  for(r in 1:r_max){
    Xmat <- rbind(Xmat, XList[[r]])
  }
  ## yangyi
  
  # require(mclust)
  message('Start computing intial values... \n')
  princ1 <- wpca(Xmat, q, weighted = wpca_int)
  hZ <- princ1$PCs
  W0 <- princ1$loadings
  set.seed(seed)
  ### parallel to obtain the initial values 
  alpha <- FALSE
  nK <- length(K)
  if(nK>1){
    ## at most use 80% CPU
    if(is.null(coreNum_int)){
      cores <- ifelse(nK < 10,  
                      nK, 10)
    }else{
      cores <- coreNum_int
    }
    
    if(Sys.info()[1]=="Windows"){
      cl <- parallel::makeCluster(cores) # memory can not be shared type in Windows.
    }else{
      cl <- parallel::makeCluster(cores, type='FORK') # memory-shared type in linux or Mac.
    }
    # Mclust <- mclust::Mclust
    # mclustBIC <- mclust::mclustBIC
    message("Starting parallel computing initial values...")
    # parallel::clusterExport(cl, varlist = c("Mclust", "mclustBIC"))
    # Run
    
    intList <- parallel::parLapply(cl, X=K, parafun_int, Z=hZ, alpha=alpha, 
                                   Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                                   int.model =int.model, verbose=verbose)
    parallel::stopCluster(cl)
  }else{
    intList <- list(parafun_int(K, Z=hZ, alpha=alpha, Sigma_equal, Sigma_diag,
                                int.model =int.model ,verbose=verbose))
  }
  
  alpha0List = list()
  Mu0List = list()
  Sigma0List = list()
  ymat = matrix(0, n, nK)
  for (kk in 1:nK){
    # print(kk)
    
    ymat[,kk] <- intList[[kk]]$ymatk
    
    alpha0List[[kk]] <- intList[[kk]]$alpha0k
    
    Mu0List[[kk]] <- intList[[kk]]$Mu0k
    
    Sigma0List[[kk]] <- intList[[kk]]$Sigma0k
  }
  
  
  
  
  message("----Fitting PRECAST model----------------\n")
  beta0  =1.5
  resList <- idrsc2Cpp(XList, AdjList, Adjlist_car, hZ, ymat, Mu0List, 
                       Sigma0List, W0, alpha0List, beta0, 
                       beta_grid, maxIter_ICM, 
                       maxIter, epsLogLik, verbose,(!error_heter),
                       Sigma_equal, Sigma_diag, mix_prop_heter, Sp2, max(K), min(K), coreNum)
  
  
  
  
  para_settings <- list(K=K, n=n, p=p,q=q,r_max=r_max,
                        Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                        mix_prop_heter=mix_prop_heter)
  attr(resList, "para_settings") <- para_settings
  class(resList) <- "SeqK_PRECAST_Object"
  return(resList)
  
  
}




idrsc <- function(XList, q, K, AdjList=NULL,  Adjlist_car=NULL, posList = NULL, platform = "ST",
                  beta_grid=seq(0.2,4, by=0.2),
                   maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE,
                   mix_prop_heter=TRUE, Sigma_equal=FALSE, Sigma_diag=TRUE,error_heter=TRUE, Sp2=TRUE,
                   wpca_int=FALSE,int.model='EEE', seed=1,coreNum = 1, coreNum_int=coreNum){
  
  
  
  
  
  
  r_max <- length(XList)
  
  pf <- unlist(lapply(XList, ncol))
  # if(is.null(posList) && is.null(AdjList)) stop("Check the values of arguments!")
  if(sum(pf!=pf[1]) > 0) stop("ncols of XList's components must be same!")
  if(is.null(posList) && is.null(AdjList)) stop("posList or AdjList must be provided one!")
  n <- sum(sapply(XList, nrow))
  p <-  pf[1]
  message("Intergrative data info.: ", r_max, ' samples, ', p , " genes X ", n, " spots------")
  message("iDR-SC model setting: ", 'error_heter=',error_heter,", Sigma_equal=",Sigma_equal,
          ", Sigma_diag=", Sigma_diag, ', mix_prop_heter=', mix_prop_heter)
  
  if(is.null(AdjList) && !is.null(posList))
    AdjList <- lapply(posList, getAdj, platform=platform)
  if(any(sapply(AdjList, nrow)!=sapply(XList, nrow)))  
    stop('The dimension of Adj matrix does not match that of X!\n')
  
  if(is.null(Adjlist_car)) Adjlist_car <- AdjList
  #XList <- lapply(XList, scale, scale=scale_flag)
  Xmat <- NULL
  for(r in 1:r_max){
    Xmat <- rbind(Xmat, XList[[r]])
  }
  
  
  # require(mclust)
  message('Start computing intial values... \n')
  princ1 <- wpca(Xmat, q, weighted = wpca_int)
  hZ <- princ1$PCs
  W0 <- princ1$loadings
  set.seed(seed)
  ### parallel to obtain the initial values 
  alpha <- FALSE
  nK <- length(K)
  if(nK>1){
    ## at most use 80% CPU
    if(is.null(coreNum_int)){
      cores <- ifelse(nK < parallel::detectCores()*0.8,  
                      nK, parallel::detectCores()*0.8)
    }else{
      cores <- coreNum_int
    }
    
    if(Sys.info()[1]=="Windows"){
      cl <- parallel::makeCluster(cores) # memory can not be shared type in Windows.
    }else{
      cl <- parallel::makeCluster(cores, type='FORK') # memory-shared type in linux or Mac.
    }
    
  message("Starting parallel computing initial values...")
  # parallel::clusterExport(cl, varlist = c("Mclust", "mclustBIC"))
  # Run
  
  intList <- parallel::parLapply(cl, X=K, parafun_int, Z=hZ, alpha=alpha, 
                                 Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                                 int.model =int.model, verbose=verbose)
  parallel::stopCluster(cl)
  }else{
    intList <- list(parafun_int(K, Z=hZ, alpha=alpha, Sigma_equal, Sigma_diag,
                                int.model =int.model ,verbose=verbose))
  }
  
  alpha0List = list()
  Mu0List = list()
  Sigma0List = list()
  ymat = matrix(0, n, nK)
  for (kk in 1:nK){
    # print(kk)
    
    ymat[,kk] <- intList[[kk]]$ymatk
    
    alpha0List[[kk]] <- intList[[kk]]$alpha0k
    
    Mu0List[[kk]] <- intList[[kk]]$Mu0k
    
    Sigma0List[[kk]] <- intList[[kk]]$Sigma0k
  }
  
  
  
    
    message("----Fitting iDR-SC model----------------\n")
    beta0  =1.5
    resList <- idrsc2Cpp(XList, AdjList, Adjlist_car, hZ, ymat, Mu0List, 
                           Sigma0List, W0, alpha0List, beta0, 
                           beta_grid, maxIter_ICM, 
                           maxIter, epsLogLik, verbose,(!error_heter),
                           Sigma_equal, Sigma_diag, mix_prop_heter, Sp2, max(K), min(K), coreNum)
      
    
    
    
    para_settings <- list(K=K, n=n, p=p,q=q,r_max=r_max,
                          Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                          mix_prop_heter=mix_prop_heter)
    attr(resList, "para_settings") <- para_settings
    class(resList) <- "SeqK_PRECAST_Object"
    return(resList)
    

}


## Define a new function mycluster to avoid using parallel::clusterExport in parallel
mycluster <- function(Z, G, int.model='EEE', verbose=FALSE){
  # require(mclust)
  mclus2 <- mclust::Mclust(Z, G=G,modelNames =int.model ,verbose=verbose)
  return(mclus2)
}
parafun_int <- function(k, Z, alpha, Sigma_equal, Sigma_diag,int.model='EEE', verbose=FALSE){
  
  mclus2 <- mycluster(Z, G=k, int.model =int.model ,verbose=verbose)
  
  ymatk <- mclus2$classification
  
  if(alpha){
    alpha0k <- mclus2$parameters$pro
  }else{
    alpha0k <- rep(0, k)
  }
  
  Mu0k <- t(mclus2$parameters$mean)
  Sigmak <- mclus2$parameters$variance$sigma
  Sigma0k <- Sigmak
  if(Sigma_diag){
    Sigma0k <- array(0, dim=dim(Sigmak))
    for(kk in 1:k){
      diag(Sigma0k[,,kk]) <- diag(Sigmak[,,kk])
    }
    Sigmak <- Sigma0k
  }
  if(Sigma_equal){
    for(kk in 1:k){
      Sigma0k[,,kk] <- apply(Sigmak, c(1,2), mean)
    }
  }
  Pi0k <- mclus2$parameters$pro
  return(list(ymatk=ymatk, alpha0k=alpha0k, Mu0k=Mu0k, Sigma0k=Sigma0k, Pi0k=Pi0k))
}


degree_freedom <- function(K, paraList){
  
  p  <- paraList$p; q <- paraList$q
  r_max <- paraList$r_max; Sigma_equal <- paraList$Sigma_equal
  Sigma_diag <- paraList$Sigma_diag;  mix_prop_heter <- paraList$mix_prop_heter
  if(!Sigma_equal){
    if(mix_prop_heter){
      
      if(Sigma_diag){ # our model
        
        # # beta_r + W + Lam_r + Mu_k + Sigma_k + Psi_r
        dfree <- 1*r_max+p*(q+r_max) + K*(q+ q) + q*(q+1)/2*r_max
      }else{
        dfree <- 1*r_max+p*(q+r_max) + K*(q + q*(q+1)/2.0)+  q*(q+1)/2*r_max
      }
    }else{
      
      if(Sigma_diag){
        
        # # beta + W + Lam_r + Mu+Sigma + tau_r
        dfree <- 1+p*(q+r_max) + K*(q+q) +  q*(q+1)/2*r_max
      }else{
        dfree <- 1+p*(q+r_max) + K*(q+q*(q+1)/2.0)+  q*(q+1)/2*r_max
      }
    }
  }else{
    if(mix_prop_heter){
      if(Sigma_diag){ ## This is our model setting
        message("Sigma is set to diagonal matrices \n")
        # # beta_r + W + Lam_r + Mu+ Sigma_r + Phi_r
        dfree <- 1*r_max+p*(q+r_max) + K*q+ q + r_max *q*(q+1)/2
      }else{
        message("Sigma is set to dense matrices \n")
        dfree <- 1*r_max+p*(q+r_max) + K*q+ q*(q+1)/2.0+ r_max *q*(q+1)/2
      }
    }else{
      if(Sigma_diag){
        message("Sigma is set to diagonal matrices \n")
        # # beta + W + Lam_r + Mu_k +Sigma_r + Phi_r
        dfree <- 1+p*(q+r_max) + K*q+ q  + r_max *q*(q+1)/2
      }else{
        message("Sigma is set to dense matrices \n")
        dfree <- 1+p*(q+r_max) + K*q+ q*(q+1)/2.0 + r_max *q*(q+1)/2
      }
    }
  }
  return(dfree)
}

selectModel <- function(obj, criteria = 'MBIC',pen_const=1, return_para_est=FALSE) UseMethod("selectModel")
selectModel.SeqK_PRECAST_Object <- function(obj, criteria = 'MBIC',pen_const=1, return_para_est=FALSE){
  if(!inherits(obj, 'SeqK_PRECAST_Object')) stop('obj must be "SeqK_PRECAST_Object"!\n')
  para_settings <- attr(obj, 'para_settings')
  K_set <- para_settings$K
  nK <- length(K_set)
  icVec <-  rep(Inf, nK)
  if(length(obj) != nK) stop("The length of obj must be equal to length of K!")
  n <- para_settings$n; p <- para_settings$p
  for(k in 1:nK){
      resList <- obj[[k]]
      dfree <- degree_freedom(K_set[k], para_settings)
      icVec[k] <- switch(criteria, 
             MAIC = -2.0* resList$loglik +dfree * 2 * log(log(p+n))*pen_const,
             AIC = -2.0* resList$loglik +dfree * 2,
             MBIC =  -2.0* resList$loglik +dfree * log(n) * log(log(p+n))*pen_const, 
             BIC =  -2.0* resList$loglik +dfree * log(n) * log(log(p+n))*pen_const )
      
    
    
  }
  bestK <- K_set[which.min(icVec)]
  icMat <- cbind(K=K_set, IC=icVec)
  resList <- obj[[which.min(icVec)]]
  cluster_PCList <- list(bestK= bestK, cluster=resList$cluster, hZ = resList$hZ, Rf = resList$Rf,
                         hV=resList$hV, hW = resList$W,icMat=icMat)
  if(return_para_est){
    attr(cluster_PCList, 'fit') <- resList[-c(1:4)]
  }
  
  return(cluster_PCList)
  
}



gendataInte_sp <- function(height1=30, width1=30,height2=height1, width2=width1, 
                           p =100, q=10, K=7,  G=4, beta=1.2, sigma2=1, 
                           alpha=8, seed=1, view=TRUE){
  # height <- 70
  # width <- 70
  # G <- 4
  # beta <- 1.0
  # K <- 7
  # q <- 10
  # p <- 1000
  if(q <2) stop("error:gendata_sp::q must be greater than 2!")
  
  if(length(beta) <2) beta <- rep(beta, 2)
  # require(GiRaF)
  # require(MASS)
  n1 <- height1 * width1 # # of cell in each indviduals 
  n2 <- height2 * width2
  index1 <- 1:n1; index2 <- (n1+1):(n1+n2)
  ## generate deterministic parameters, fixed after generation
 
  if(length(sigma2)==1){
    Lambda1 <- sigma2*(abs(rnorm(p, sd=1)))
    Lambda2 <- sigma2*(abs(rnorm(p, sd=1)))
  }else{
    Lambda1 <- rep(sigma2[1], p)
    Lambda2 <- rep(sigma2[2], p)
  }
  
  
  W1 <- matrix(rnorm(p*q), p, q)
  W <- qr.Q(qr(W1))
  
  mu <- matrix(0, q,  K)
  diagmat = array(0, dim = c(q, q, K))
  if(q > K){
    q1 <- floor(K/2)
    for(j in 1:q1){
      if(j <= (q1/2)) mu[j,j] <- alpha
      if(j > (q1/2)) mu[j,j] <- -alpha
    }
    mu[(q1+1):q, K] <- -alpha
    
  }else if(q <= K){
    for(k in 1:K)
      mu[,k] <- rep(alpha/8 *k, q) #
  }
  for(k in 1:K){
    tmp  <- rep(1, q)
    if(k <= K/2){
      tmp[q] <- alpha
    }
    diag(diagmat[,,k]) <- tmp
  }
  
  Mu <- t(mu)
  Sigma <- diagmat
  tau2 <- rep(1, q); tau2[1] <- 4; tau2[2] <- -4
  set.seed(seed)
  # generate the spatial dependence for state variable y, a hidden Markov RF
  y1 <- sampler.mrf(iter = n1, sampler = "Gibbs", h = height1, w = width1, ncolors = K, nei = G, param = beta[1],
                    initialise = FALSE, view = view)
  y1 <- c(y1) + 1
  y2 <- sampler.mrf(iter = n2, sampler = "Gibbs", h = height2, w = width2, ncolors = K, nei = G, param = beta[2],
                    initialise = FALSE, view = view)
  y2 <- c(y2) + 1
  
  Z1 <- matrix(0, n1, q)
  Z2 <- matrix(0, n2, q)
  for(k in 1:K){
    nk1 <- sum(y1==k)
    Z1[y1==k, ] <- MASS::mvrnorm(nk1, Mu[k,], Sigma[,,k])
    
    nk2 <- sum(y2==k)
    Z2[y2==k, ] <- MASS::mvrnorm(nk2, Mu[k,]+tau2, Sigma[,,k])
  }
  X1 <- Z1 %*% t(W)+  MASS::mvrnorm(n1, rep(0,p), diag(Lambda1))
  X2 <- Z2 %*% t(W) + MASS::mvrnorm(n2, rep(0,p), diag(Lambda2))
  # compute the strength of signal
  svd_Sig1 <- svd(cov(Z1))
  W11 <- W %*% svd_Sig1$u %*% diag(sqrt(svd_Sig1$d))
  snr1 <- sum(svd(W11)$d^2) / (sum(svd(W11)$d^2)+ sum(Lambda1))
  svd_Sig2 <- svd(cov(Z2))
  W22 <- W %*% svd_Sig2$u %*% diag(sqrt(svd_Sig2$d))
  snr2 <- sum(svd(W22)$d^2) / (sum(svd(W22)$d^2)+ sum(Lambda2))
  
  message("SNR1=", round(snr1,4),", SNR2=", round(snr2, 4),  '\n')
  
  # make position
  pos1 <- cbind(rep(1:height1, width1), rep(1:height1, each=width1))
  pos2 <- cbind(rep(1:height2, width2), rep(1:height2, each=width2))
  return(list(X1=X1,X2=X2, Z1=Z1, Z2=Z2, cluster=c(y1,y2), Mu=Mu, Sigma=Sigma,
              W=W, Lambda1=Lambda1, Lambda2=Lambda2, tau2=tau2, 
              beta=beta,pos1=pos1,pos2=pos2, snr=c(snr1,snr2),
              index1=index1,index2=index2))
  
}



find_neighbors <- function(pos, platform=c('ST', "Visium")) {
  # require(purrr)
  # require(S4Vectors)
  if (tolower(platform) == "visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (tolower(platform) == "st") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  colnames(pos) <- c("row", "col")
  pos <- DataFrame(pos)
  spot.positions <- pos
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  #neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  ## df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  df_j <- lapply(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          nrow(pos), " spots.")
  
  n <- length(df_j) 
  
  D <- matrix(0,  nrow = n, ncol = n)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  ij <- which(D != 0, arr.ind = T)
  ij
}
getAdj_reg <- function(pos, platform= "Visium"){
  ij <- find_neighbors(pos, platform)
  # library(Matrix)
  n <- nrow(pos)
  Adj_sp <- Matrix::sparseMatrix(ij[,1], ij[,2], x = 1, dims=c(n, n))
  return(Adj_sp)
}
# getAdj <- function(pos, platform="Visisum"){
#   
#   if(tolower(platform) %in% tolower(c('ST', "Visium"))){
#     Adj_sp <- getAdj_reg(pos, platform)
#   }else{
#     Adj_sp <- getAdj_auto(pos)
#   }
# }

# getAdj_auto <- function(pos){
#   if (!inherits(pos, "matrix"))
#     stop("method is only for  matrix object!")
#   
#   radius.lower <- 1
#   radius.upper <- 50
#   start.radius <- 1
#   Med <- 0
#   while(!(Med >= 4 && Med <=6)){
#     
#     Adj_sp <- getneighborhood_fast(pos, radius=start.radius)
#     Med <- summary(Matrix::rowSums(Adj_sp))['Median']
#     if(Med < 4){
#       radius.lower <- start.radius
#       start.radius <- (radius.lower + radius.upper)/2
#     }else if(Med >6){
#       radius.upper <- start.radius
#       start.radius <- (radius.lower + radius.upper)/2
#     }
#     message("Current radius is ", start.radius, '\n')
#     message("Median of neighbors is ", Med, '\n')
#   }
#   return(Adj_sp)
# }
# 

# getAdj_manual <- function(pos, radius){
#   if (!inherits(pos, "matrix"))
#     stop("pos must be a matrix!")
#   if(radius <=0){
#     stop('radius must be a positive real!')
#   }
#   
#   Adj_sp <- getneighborhood_fast(pos, radius=radius)
#   return(Adj_sp)
# }

getAdj_fixedNumber <- function(pos, number=6){
  if(nrow(pos)< 6*4) stop("nrow of pos must be greater than 4*number!")
  # require(Matrix)
  Adj_sp <- get_fixedNumber_neighbors(pos, number)
  return(Adj_sp)
}

wpca <- function(X, q, weighted=T){
  if(!is.matrix(X)) stop("wpca: X must be a matrix!")
  if(q< 1) stop("wpca: q must be a positive integer!")
  X <- scale(X, scale=F) # centralize
  out <- wpcaCpp(X, q, weighted)
  return(out)
}





# Used in real data analysis ----------------------------------------------



## Get basic information
## format the log-normalized gene expression based on the selected genes.
getXList <- function(seuList, genelist){
  # require(Seurat)
  
  seuList <- lapply(seuList, NormalizeData)
  XList <- list()
  nsample <- length(seuList)
  indexList <- list()
  nr <- 0
  posList <- list()
  for(i in 1:nsample){
    XList[[i]] <- Matrix::t(seuList[[i]]@assays$RNA@data[genelist,])
    indexList[[i]] <- (nr+1):(nrow(XList[[i]] )+nr)
    nr <- nr + nrow(XList[[i]] )
    posList[[i]] <- cbind(seuList[[i]]$row, seuList[[i]]$col)
  }
  
  return(list(XList=XList, posList=posList, 
              indxList=indexList))
}



