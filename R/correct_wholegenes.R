### These utility functions are used to correct the whole genes
mat2list <- function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

vec2list <- function(y_int, nvec){
  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}

get_correct_exp <- function(XList, RfList,  houseKeep, covariateList=NULL, q_unwanted=10, subsample_rate=NULL, sample_seed=1){
  
  
  if(!all(sapply(XList, is.matrix))){
    XList <- lapply(XList, as.matrix)
  }
  nvec <- sapply(XList, nrow)
  XList_sub <- pbapply::pblapply(XList, function(x) x[,houseKeep])
  M0 <- wpca(matlist2mat(XList_sub), q=q_unwanted, FALSE)$PCs 
  HList <- mat2list(M0, nvec=nvec)
  
  Rf <- matlist2mat(RfList)
  colnames(Rf) <- paste0("Rf",1:ncol(Rf))
  if(!is.null(covariateList)){
    # covariates <- matlist2mat(covariateList)
    # covariates <- as.matrix(covariates)
    covarites_df <- dfList2df(covariateList)
    covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
    rm(covariateList, covarites_df)
    Rf <- cbind(Rf, covariates[,-1])
    rm(covariates)
  }
  rm(RfList)
  
  # subsampling to speed up the computation!
  if(is.null(subsample_rate)) subsample_rate <- 1
  index_List <- get_indexList(XList)
  set.seed(sample_seed)
  index_subsample <- sort(sample(sum(nvec), floor(sum(nvec)*subsample_rate)))
  ## calculate the number of indices belonging to the index of each slide
  nvec_subsample <- rep(NA, length(nvec))
  for(i in 1:length(nvec_subsample)){
    ## message("i = ", i)
    nvec_subsample[i] <- sum(index_subsample%in% index_List[[i]])
  }
  index_List_new <- lapply(XList, function(x) 1: nrow(x))
  index_subsample_new <- unlist(index_List_new)[index_subsample]
  index_subsampleList <- vec2list(index_subsample_new, nvec_subsample)
  RfList <- mat2list(Rf, nvec=nvec)
  
  XList_sub <- list(); RList_sub <- list(); HList_sub <- list()
  AdjList_sub <- list()
  for(i in 1:length(XList)){
    # message("i = ", i)
    index_tmp <- index_subsampleList[[i]]
    XList_sub[[i]] <- XList[[i]][index_tmp,]
    RList_sub[[i]] <- RfList[[i]][index_tmp, ]
    HList_sub[[i]] <- HList[[i]][index_tmp, ]
  }
  rm(RfList,  HList, index_tmp)
  
  
  ### XList <-  lapply(XList, scale, scale=FALSE)
  X_sub <- matlist2mat(XList_sub)
  rm(XList_sub)
  Rf <- matlist2mat(RList_sub)
  colnames(Rf) <- NULL
  H <- matlist2mat(HList_sub)
  colnames(H) <- NULL
  nc_M0 <- ncol(H)
  lm1 <- lm(X_sub~ 0+ cbind(H, Rf))
  coefmat <- coef(lm1)[c(1:nc_M0),]
  #row.names(coef(lm1))
  rm(lm1, X_sub)
  hX <- matlist2mat(XList) - M0 %*% coefmat
  return(hX)
}

get_correct_mean_exp <- function(XList,  hVList, covariateList=NULL, subsample_rate=NULL, sample_seed=1){
  
  ## XList <- lapply(XList, scale, scale=FALSE)
  
  if(!all(sapply(XList, is.matrix))){
    XList <- lapply(XList, as.matrix)
  }
  
  r_max <- length(XList)
  X0 <- XList[[1]]
  hV0 <- hVList[[1]]
  if(r_max>1){
    for(r in 2:r_max){
      X0 <- rbind(X0, XList[[r]])
      hV0 <- rbind(hV0, hVList[[r]])
      
    }
  }
  
  rm(XList)
  if(!is.null(covariateList)){
    # covariates <- matlist2mat(covariateList)
    # covariates <- as.matrix(covariates)
    # covariates <- cbind(1, covariates)
    covarites_df <- dfList2df(covariateList)
    covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
    rm(covariateList, covarites_df)
  }else{
    covariates <- matrix(1, nrow=nrow(hV0), ncol=1)
  }
  nc_M0 <- ncol(hV0)
  if(!is.null(subsample_rate)){
    set.seed(sample_seed)
    id_sub <- sort(sample(nrow(X0), floor(nrow(X0)*subsample_rate)))
    X_sub <- X0[id_sub,]; hV_sub <- hV0[id_sub,]; covariates_sub <- covariates[id_sub,]
    lm1 <- lm(X_sub~  0+ cbind(hV_sub, covariates_sub))
  }else{
    lm1 <- lm(X0~  0+ cbind(hV0, covariates))
  }
  rm(covariates)
  coefmat <- coef(lm1)[c(1:nc_M0),]
  # row.names(coef(lm1))
  rm(lm1)
  X0 - hV0 %*% coefmat
  
  
}


IntegrateSpaData <- function(PRECASTObj, species="Human", custom_housekeep=NULL,
                             covariates_use=NULL, seuList=NULL, subsample_rate=1, sample_seed=1){
  # suppressMessages(require(Matrix))
  # suppressMessages(require(Seurat))
  
  verbose <- PRECASTObj@parameterList$verbose
  if(!inherits(PRECASTObj, "PRECASTObj")) 
    stop("IntegrateSpaData: Check the argument: PRECASTObj!  PRECASTObj must be a PRECASTObj object.")
  if(is.null(PRECASTObj@seulist)) 
    stop("IntegrateSpaData: Check the argument: PRECASTObj! The slot seulist in PRECASTObj is NULL!")
  
  if(length(subsample_rate)>1 | subsample_rate<0 | subsample_rate>1) stop("subsample_rate must be a real between 0 and 1")
  
  if(!tolower(species) %in% c("human", "mouse", "unknown")) 
    stop("IntegrateSpaData: Check the argument: species! it must be one of 'Human', 'Mouse' and 'Unknown'!")
  
  defAssay_vec <- sapply(PRECASTObj@seulist, DefaultAssay)
  if(any(defAssay_vec!=defAssay_vec[1])) warning("IntegrateSpaData: there are different default assays in PRECASTObj@seulist that will be used to integrating!")
  n_r <- length(defAssay_vec)
  if(is.null(seuList) && (!is.null(PRECASTObj@seuList))){
    message("Use PRECASTObj@seuList as input to remove unwanted variation since it is preserved in PRECASTObj.")
    seuList <- PRECASTObj@seuList
  }
  if(!is.null(seuList)){
    message("seuList is not NULL. Filter the spots in seuList but not in PRECASTObj!")
    if(!is.null(names(seuList)) && !is.null(names(PRECASTObj@seulist))){
      if(!all(names(seuList) == names(PRECASTObj@seulist))){
        stop("The names of seuList must be the same as the names of PRECASTObj@seulist.")
      }
    }
    seuList <- lapply(seq_along(PRECASTObj@seulist), function(j){
      # j <- 1
      seu1 <- PRECASTObj@seulist[[j]]
      seu <- seuList[[j]]
      if(length(setdiff(colnames(seu1), colnames(seu)))){
        stop("The spot's name in PRECASTObj@seulist[[",j, "]] must be a subset of the spot's name in that of seuList!")
      }
      seu <- seu[, colnames(seu1)]
      return(seu)
    })
    seuList <- lapply(seuList, NormalizeData, verbose=FALSE)
    XList <- lapply(1:n_r,  function(r) Matrix::t(get_data_fromSeurat(seuList[[r]], assay = defAssay_vec[r], slot='data')))
  }else{
    XList <- lapply(1:n_r,  function(r) Matrix::t(get_data_fromSeurat(PRECASTObj@seulist[[r]], assay = defAssay_vec[r], slot='data')))
    
  }
  
  
  
  if(!is.null(covariates_use)){
    # covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[, covariates_use])
    covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[covariates_use])
  }else{
    covariateList <- NULL
  }
  
  if(tolower(species) =='mouse'){
    for(r in 1:n_r){
      colnames(XList[[r]]) <- firstup(colnames(XList[[r]]))
    }
    if(!is.null(custom_housekeep)){
      custom_housekeep <- firstup(custom_housekeep)
    }
  }
  if(tolower(species) =='human'){
    for(r in 1:n_r){
      colnames(XList[[r]]) <- toupper(colnames(XList[[r]]))
    }
    if(!is.null(custom_housekeep)){
      custom_housekeep <- toupper(custom_housekeep)
    }
  }
  
  barcodes_all <- lapply(XList, row.names)
  if(any(duplicated(unlist(barcodes_all)))){
    
    for(r in 1:n_r){
      row.names(XList[[r]]) <- paste0(row.names(XList[[r]]), r)
    }
  }
  genelist <- colnames(XList[[1]])
  lower_species <- tolower(species) 
  houseKeep <- switch (lower_species,
                       human = {
                         # data(Human_HK_genes)
                         intersect(toupper(genelist), PRECAST::Human_HK_genes$Gene)
                       },
                       mouse={
                         #data(Mouse_HK_genes)
                         intersect(firstup(genelist), PRECAST::Mouse_HK_genes$Gene)
                       },
                       unknown={
                         character()
                       }
  )
  houseKeep <- c(houseKeep, custom_housekeep)
  houseKeep <- intersect(houseKeep, colnames(XList[[1]]))
  
  nvec <- sapply(XList, nrow)
  if(sum(nvec)> 8e4){
    
    subsample_rate <- 5e4/sum(nvec)
    message("The total number of spots exceeds 8e4, thus the subsampling schema will be used to speed up computation.")
  }
  if(subsample_rate <1 && subsample_rate>0) message("IntegrateSRTData: the subsampling schema will be used to speed up computation since subsample_rate is smaller than 1.")
  
  
  tstart <- Sys.time()
  if(length(houseKeep) < 5){
    if(verbose){
      message("Using only PRECAST results to obtain the batch corrected gene expressions since species is unknown or the genelist in PRECASTObj has less than 5 overlapp with the housekeeping genes of given species.")
      message("Start integration...")
    }
    
    hX <- get_correct_mean_exp(XList, PRECASTObj@resList$hV,  covariateList=covariateList,
                                subsample_rate=subsample_rate,sample_seed=sample_seed)
  }else{
    if(verbose){
      message("Using bouth housekeeping gene and PRECAST results to obtain the batch corrected gene expressions.")
      message("Start integration...")
    }
    # matlist2mat <- PRECAST:::matlist2mat
    # get_indexList <- PRECAST:::get_indexList
    hX <- get_correct_exp(XList, PRECASTObj@resList$Rf, houseKeep=houseKeep, q_unwanted=min(10, length(houseKeep)),
                           covariateList=covariateList,subsample_rate=subsample_rate,sample_seed=sample_seed)
    # hX[1:4,1:4]
  }
  .logDiffTime(sprintf(paste0("%s Data integration finished!"), "*****"), t1 = tstart, verbose = verbose)
  
  if(verbose)
    message("Put the data into a new Seurat object...")
  tstart <- Sys.time()
  meta_data <- data.frame(batch=factor(get_sampleID(XList)), cluster= factor(unlist(PRECASTObj@resList$cluster)))
  row.names(meta_data) <- row.names(hX)
  rm(XList)
  count <- sparseMatrix(i=1,j=1, x=0, dims=dim(t(hX)))
  row.names(count) <- colnames(hX)
  colnames(count) <- row.names(hX)
  seuInt <- CreateSeuratObject(counts = count, assay = 'PRE_CAST', meta.data=meta_data)
  if(inherits(seuInt[['PRE_CAST']], "Assay5") ){
    
    seuInt <- SetAssayData(object = seuInt, layer='data', assay = "PRE_CAST", new.data =  t(hX))
  }else{
    seuInt[['PRE_CAST']]@data <- t(hX)
  }
  
  rm(hX)
  
  # seuInt <- CreateSeuratObject(assay, meta.data=meta_data, assay = 'PRECAST')
  seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$hZ), seuInt, embed_name = 'PRECAST', assay='PRE_CAST')
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  seuInt <- Add_embed(matlist2mat(posList), seuInt, embed_name = 'position', assay='PRE_CAST')
  Idents(seuInt) <- factor(meta_data$cluster)
  
  .logDiffTime(sprintf(paste0("%s New Seurat object is generated!"), "*****"), t1 = tstart, verbose = verbose)
  
  
  return(seuInt)
}


# 
# get_correct_exp <- function(XList, RfList,  houseKeep, covariateList=NULL, q_unwanted=10){
#   
#   
#   if(!all(sapply(XList, is.matrix))){
#     XList <- lapply(XList, as.matrix)
#   }
#   XList_sub <- pbapply::pblapply(XList, function(x) x[,houseKeep])
#   M0 <- wpca(matlist2mat(XList_sub), q=q_unwanted, FALSE)$PCs 
#   
#   
#   Rf <- matlist2mat(RfList)
#   colnames(Rf) <- paste0("Rf",1:ncol(Rf))
#   rm(RfList)
#   if(!is.null(covariateList)){
#     # covariates <- matlist2mat(covariateList)
#     # covariates <- as.matrix(covariates)
#     covarites_df <- dfList2df(covariateList)
#     covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
#     rm(covariateList, covarites_df)
#     Rf <- cbind(Rf, covariates[,-1])
#     rm(covariates)
#   }
#   ### XList <-  lapply(XList, scale, scale=FALSE)
#   X0 <- matlist2mat(XList)
#   rm(XList)
#   nc_M0 <- ncol(M0)
#   lm1 <- lm(X0~ 0+ cbind(M0, Rf))
#   coefmat <- coef(lm1)[c(1:nc_M0),]
#   #row.names(coef(lm1))
#   rm(lm1)
#   hX <- X0 - M0 %*% coefmat
#   return(hX)
# }
# 
# get_correct_mean_exp <- function(XList,  hVList, covariateList=NULL){
#   
#   ## XList <- lapply(XList, scale, scale=FALSE)
#   
#   if(!all(sapply(XList, is.matrix))){
#     XList <- lapply(XList, as.matrix)
#   }
#   
#   r_max <- length(XList)
#   X0 <- XList[[1]]
#   hV0 <- hVList[[1]]
#   if(r_max>1){
#     for(r in 2:r_max){
#       X0 <- rbind(X0, XList[[r]])
#       hV0 <- rbind(hV0, hVList[[r]])
#       
#     }
#   }
#   
#   rm(XList)
#   if(!is.null(covariateList)){
#     # covariates <- matlist2mat(covariateList)
#     # covariates <- as.matrix(covariates)
#     # covariates <- cbind(1, covariates)
#     covarites_df <- dfList2df(covariateList)
#     covariates <-  model.matrix.lm(object = ~.+1, data = covarites_df, na.action = "na.pass")
#     rm(covariateList, covarites_df)
#   }else{
#     covariates <- matrix(1, nrow=nrow(hV0), ncol=1)
#   }
#   
#   nc_M0 <- ncol(hV0)
#   lm1 <- lm(X0~  0+ cbind(hV0, covariates))
#   coefmat <- coef(lm1)[c(1:nc_M0),]
#   # row.names(coef(lm1))
#   rm(lm1)
#   X0 - hV0 %*% coefmat
#   
#   
# }
# 
# 
# IntegrateSpaData <- function(PRECASTObj, species="Human", custom_housekeep=NULL,
#                              covariates_use=NULL){
#   # suppressMessages(require(Matrix))
#   # suppressMessages(require(Seurat))
#   
#   verbose <- PRECASTObj@parameterList$verbose
#   if(!inherits(PRECASTObj, "PRECASTObj")) 
#     stop("IntegrateSpaData: Check the argument: PRECASTObj!  PRECASTObj must be a PRECASTObj object.")
#   if(is.null(PRECASTObj@seulist)) 
#     stop("IntegrateSpaData: Check the argument: PRECASTObj! The slot seulist in PRECASTObj is NULL!")
#   
#   
#   if(!tolower(species) %in% c("human", "mouse", "unknown")) 
#     stop("IntegrateSpaData: Check the argument: species! it must be one of 'Human', 'Mouse' and 'Unknown'!")
#   
#   defAssay_vec <- sapply(PRECASTObj@seulist, DefaultAssay)
#   if(any(defAssay_vec!=defAssay_vec[1])) warning("IntegrateSpaData: there are different default assays in PRECASTObj@seulist that will be used to integrating!")
#   n_r <- length(defAssay_vec)
#   
#   XList <- lapply(1:n_r,  function(r) Matrix::t(GetAssayData(PRECASTObj@seulist[[r]], assay = defAssay_vec[r], slot='data')))
#   
#   if(!is.null(covariates_use)){
#     # covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[, covariates_use])
#     covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[covariates_use])
#   }else{
#     covariateList <- NULL
#   }
#   
#   if(tolower(species) =='mouse'){
#     for(r in 1:n_r){
#       colnames(XList[[r]]) <- firstup(colnames(XList[[r]]))
#     }
#     if(!is.null(custom_housekeep)){
#       custom_housekeep <- firstup(custom_housekeep)
#     }
#   }
#   if(tolower(species) =='human'){
#     for(r in 1:n_r){
#       colnames(XList[[r]]) <- toupper(colnames(XList[[r]]))
#     }
#     if(!is.null(custom_housekeep)){
#       custom_housekeep <- toupper(custom_housekeep)
#     }
#   }
#   
#   barcodes_all <- lapply(XList, row.names)
#   if(any(duplicated(unlist(barcodes_all)))){
#     
#     for(r in 1:n_r){
#       row.names(XList[[r]]) <- paste0(row.names(XList[[r]]), r)
#     }
#   }
#   genelist <- colnames(XList[[1]])
#   lower_species <- tolower(species) 
#   houseKeep <- switch (lower_species,
#                        human = {
#                          # data(Human_HK_genes)
#                          intersect(toupper(genelist), PRECAST::Human_HK_genes$Gene)
#                        },
#                        mouse={
#                          #data(Mouse_HK_genes)
#                          intersect(firstup(genelist), PRECAST::Mouse_HK_genes$Gene)
#                        },
#                        unknown={
#                          character()
#                        }
#   )
#   houseKeep <- c(houseKeep, custom_housekeep)
#   houseKeep <- intersect(houseKeep, colnames(XList[[1]]))
#   
#   tstart <- Sys.time()
#   if(length(houseKeep) < 5){
#     if(verbose){
#       message("Using only PRECAST results to obtain the batch corrected gene expressions since species is unknown or the genelist in PRECASTObj has less than 5 overlapp with the housekeeping genes of given species.")
#       message("Start integration...")
#     }
#     
#     hX <- get_correct_mean_exp(XList,PRECASTObj@resList$hV,  covariateList=covariateList)
#   }else{
#     if(verbose){
#       message("Using bouth housekeeping gene and PRECAST results to obtain the batch corrected gene expressions.")
#       message("Start integration...")
#     }
#     
#     hX <- get_correct_exp(XList, PRECASTObj@resList$Rf, houseKeep=houseKeep, q_unwanted=min(10, length(houseKeep)), covariateList=covariateList)
#   }
#   .logDiffTime(sprintf(paste0("%s Data integration finished!"), "*****"), t1 = tstart, verbose = verbose)
#   
#   if(verbose)
#     message("Put the data into a new Seurat object...")
#   tstart <- Sys.time()
#   meta_data <- data.frame(batch=factor(get_sampleID(XList)), cluster= factor(unlist(PRECASTObj@resList$cluster)))
#   row.names(meta_data) <- row.names(hX)
#   rm(XList)
#   count <- sparseMatrix(i=1,j=1, x=0, dims=dim(t(hX)))
#   row.names(count) <- colnames(hX)
#   colnames(count) <- row.names(hX)
#   seuInt <- CreateSeuratObject(counts = count, assay = 'PRE_CAST', meta.data=meta_data)
#   if(inherits(seuInt[['PRE_CAST']], "Assay5") ){
#     
#     seuInt <- SetAssayData(object = seuInt, slot='data', assay = "PRE_CAST", new.data =  t(hX))
#   }else{
#     seuInt[['PRE_CAST']]@data <- t(hX)
#   }
#   
#   rm(hX)
#   
#   # seuInt <- CreateSeuratObject(assay, meta.data=meta_data, assay = 'PRECAST')
#   seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$hZ), seuInt, embed_name = 'PRECAST', assay='PRE_CAST')
#   posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
#   seuInt <- Add_embed(matlist2mat(posList), seuInt, embed_name = 'position', assay='PRE_CAST')
#   Idents(seuInt) <- factor(meta_data$cluster)
#   
#   .logDiffTime(sprintf(paste0("%s New Seurat object is generated!"), "*****"), t1 = tstart, verbose = verbose)
#   
#   
#   return(seuInt)
# }
