
addnames <- function(XList){
  
  for(j in 1:length(XList)){
    row.names(XList[[j]]) <- paste0("cell",j,"_", 1:nrow(XList[[j]]))
    colnames(XList[[j]]) <- paste0("gene", 1:ncol(XList[[j]]))
    
  }
  return(XList)
}

Add_embed <- function(embed, seu, embed_name='tSNE' , assay = "RNA"){
  row.names(embed) <- colnames(seu)
  colnames(embed) <- paste0(embed_name, 1:ncol(embed))
  seu@reductions[[embed_name]] <- CreateDimReducObject(embeddings = embed, 
                                                          key = paste0(toupper(embed_name),"_"), assay = assay)
  seu
}

AddUMAP <- function(seuInt, n_comp=3, reduction='PRECAST', assay='PRE_CAST', seed=1){
  if(is.null(seuInt@reductions[[reduction]])) 
    stop(paste0("The  reduction  ", reduction, " does not exist, please change reduction or run IntegrateSpaData first!") )
  set.seed(seed)
  hZ_umap <- scater::calculateTSNE(t(seuInt@reductions[[reduction]]@cell.embeddings), ncomponents=n_comp)
  if(n_comp==3){
    seuInt <- Add_embed(hZ_umap, seuInt, embed_name = "UMAP3", assay = assay)
  }else if(n_comp==2){
    seuInt <- Add_embed(hZ_umap, seuInt, embed_name = "UMAP", assay = assay)
  }
  return(seuInt)
}

AddTSNE <- function(seuInt, n_comp=3, reduction='PRECAST', assay='PRE_CAST', seed=1){
  if(is.null(seuInt@reductions[[reduction]])) stop("The  reduction does not exist, please run IntegrateSpaData first!")
  set.seed(seed)
  hZ_tsne <- scater::calculateTSNE(t(seuInt@reductions[[reduction]]@cell.embeddings), ncomponents=n_comp)
  if(n_comp==3){
    seuInt <- Add_embed(hZ_tsne, seuInt, embed_name = "tSNE3", assay = assay)
  }else if(n_comp==2){
    seuInt <- Add_embed(hZ_tsne, seuInt, embed_name = "tSNE", assay = assay)
  }
  return(seuInt)
}
replaceStr <- function(strVec, pattern="_tmp", by="RNA"){
  n_str <- length(strVec)
  ystr <- strVec
  for(i in 1:n_str){
    ystr[i] <- sub(pattern, by, strVec[i])
  }
  return(ystr)
}
gendata_seulist <- function(height1=30, width1=30,height2=height1, width2=width1, 
                           p =100, q=10, K=7,  G=4, beta=1.2, sigma2=1, 
                           alpha=8, seed=1, view=FALSE){
  
  #height1=30; width1=30;height2=height1; width2=width1; 
  #p =100; q=10; K=7;  G=4; beta=1.2; sigma2=1; alpha=8; seed=1; view=TRUE
  
  # suppressMessages(require(GiRaF))
  # suppressMessages(require(MASS))
  # suppressMessages(require(Seurat))
  
  
  
  ###############Model parameters setting
  if(q <2) stop("error:gendata_sp::q must be greater than 1!")
  if(length(beta) <2) beta <- rep(beta, 2)
  n1 <- height1 * width1 # # of cell in each indviduals 
  n2 <- height2 * width2
  n_vec <- c(n1, n2)
  sigmaW=c(0.5,0.8,1);
  sigmaZ = 2*c(1,2, 0.5);
  qvec=rep(2, 3); 
  ## generate deterministic parameters, fixed after generation
  
  if(length(sigma2)==1){
    Lambda1 <- sigma2*(abs(rnorm(p, sd=1)))
    Lambda2 <- sigma2*(abs(rnorm(p, sd=1)))
  }else{
    Lambda1 <- rep(sigma2[1], p)
    Lambda2 <- rep(sigma2[2], p)
  }
  LamMat <- rbind(Lambda1, Lambda2)
  
  W1 <- matrix(rnorm(p*q), p, q)
  W <- qr.Q(qr(W1))
  Wlist <- list()
  for(r in 1:2){
     # sigma12 control the correlation strength, if sigma12 is small, correlation is strong
    Wtt1 <- matrix(rnorm(p* qvec[r]), p, qvec[r]) 
    W1 <- Wtt1 + sigmaW[r] * matrix(rnorm(p* qvec[r]), p, qvec[r])
    W1 <- qr.Q(qr(W1))
    Wlist[[r]] <- W1
    # cat('cor(W,W1)=', mean(cancor(W, Wlist[[r]])$cor), '\n')
  }
  
  
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
  
  
  tauMat <- matrix(0, 2, q)
  tauMat[1, ] <- rep(5, q);
  tauMat[2,1] <- 10; tauMat[2,2] <- -10 

  
  tau0Mat <- matrix(NA, 2, p)
  for(r in 1:2){
    
    tau0 <- rnorm(p, sd=2)
    tau0Mat[r, ] <- tau0
  }
  
  
  
  
  ################## Start to generate data
  set.seed(seed)
  # generate the spatial dependence for state variable y, a hidden Markov RF
  y1 <- sampler.mrf(iter = n1, sampler = "Gibbs", h = height1, w = width1, ncolors = K, nei = G, param = beta[1],
                    initialise = FALSE, view = view)
  y1 <- c(y1) + 1
  y2 <- sampler.mrf(iter = n2, sampler = "Gibbs", h = height2, w = width2, ncolors = K, 
                    nei = G, param = beta[2],
                    initialise = FALSE, view = view)
  y2 <- c(y2) + 1
  
  yList <- list(y1, y2)
  # make position
  pos1 <- cbind(rep(1:height1, width1), rep(1:height1, each=width1))
  pos2 <- cbind(rep(1:height2, width2), rep(1:height2, each=width2))
  posList <- list(pos1, pos2)
  
  Zlist <- list()
  VList <- list()
  for(r in 1:2){
    Z_tmp <- matrix(0, n_vec[r], q)
    for(k in 1:K){
      nk <- sum(yList[[r]]==k) # include conditional and sequencing batch
      if(nk > 0)
        Z_tmp[yList[[r]]==k, ] <- MASS::mvrnorm(nk, Mu[k,],Sigma[,,k])
    }
    Zlist[[r]] <- Z_tmp
    VList[[r]] <- MASS::mvrnorm(n_vec[r], tauMat[r, ],Sigma[,,1]+ 3*(r-1)*diag(q))
  }
  
  
  ## batch effect
  Zrlist <- list()
  for(r in 1:2){
    
    Zrlist[[r]] <- matrix(rnorm(n_vec[r]* qvec[r], sd=sigmaZ[r]), n_vec[r], qvec[r])
  }
  sapply(Zrlist, dim)
  
  
  
  XtList <- list()
  for(r in 1:2){
    X1 <- (Zlist[[r]] + VList[[r]] ) %*% t(W) + Zrlist[[r]] %*% t(Wlist[[r]])+ 
      MASS::mvrnorm(n_vec[r], rep(0,p), diag(LamMat[r,]))
    
    tauMat0 <- matrix(tau0Mat[r, ], n_vec[r], p, byrow = T)
    Eta <- exp((X1 + tauMat0))
    summary(colSums(Eta))
    #X1 <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
    XtList[[r]] <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
  }
  
  
  # XList <- lapply(XtList, function(x) log(1+x))
  
 
  
  
  
  XtList <- addnames(XtList)
  seulist <- list()
  for(r in 1:length(XtList)){
    
    seulist[[r]] <- CreateSeuratObject(counts= t(XtList[[r]]), assay='RNA')
    seulist[[r]]$row <- posList[[r]][,1]
    seulist[[r]]$col <- posList[[r]][,2]
    seulist[[r]]$true_cluster <- yList[[r]]
    seulist[[r]] <- Add_embed(Zlist[[r]],seulist[[r]],embed_name='trueEmbed' , assay = "RNA" )
    
    #colnames(seulist[[r]]@meta.data)[c(2,3)] <- replaceStr(colnames(seulist[[r]]@meta.data)[c(2,3)] ,"tmp", by="RNA")
  }
  # paraList <- list(LamList=list(Lambda1, Lambda2), betaVec=beta)
  # attr(seulist,"paraList") <- paraList
  
  
  return(seulist)
  
}




filter_spot <- function(seu, min_feature=0, assay="RNA"){ # each spots at least include 1 non-zero features
  col_name <- paste0("nFeature_",assay)
  idx <- seu@meta.data[,col_name] > min_feature
  seu[, idx]
  # subset(seu, subset = nFeature_RNA > min_feature)
}
# filter_spot(seu, assay='PRE_CAST')
filter_gene <- function(seu, min_spots=20, assay="RNA"){
  if(sum(dim(seu[[assay]]@counts))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@counts>0)>min_spots
  }else if(sum(dim(seu[[assay]]@data))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@data>0)>min_spots
  }else{
    stop("filter_gene: Seuat object must provide slots count or data in assay!")
  }
  
  seu[gene_flag, ]
}

## select the features for multiple samples based on a rank rule.
selectIntFeatures <- function(seulist, spaFeatureList, IntFeatures=2000){
  ## This function is used for selecting common informative features
  if(any(sapply(spaFeatureList, length)< IntFeatures))
    stop("Feature list exists number of features less than IntFeatures!")
  geneUnion <- unique(unlist(lapply(spaFeatureList, function(x) x[1:IntFeatures])))
  ## ensure each seuobject has the genes in geneUnion
  gene_delete <- unique(unlist(lapply(seulist, function(x) setdiff(geneUnion, row.names(x)))))
  geneUnion <- setdiff(geneUnion, gene_delete)
  
  
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(lapply(seulist, function(x) 
    geneUnion[pbapply::pbapply(x@assays$RNA@counts[geneUnion,],1, sd)==0])))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  
  # sort by number of datasets that identified this gene as SVG.
  nsample <- length(seulist)
  numVec <- rep(0, length(gene_Var))
  rankMat <-matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], spaFeatureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(spaFeatureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
    
  }
  
  cutNum <- sort(numVec, decreasing = T)[min(IntFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(IntFeatures, length(numVec)) - length(genelist1)
  
  gene2 <- gene_Var[numVec==cutNum]
  ### select top 2000 genes that rank 
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)
  
  return(genelist)
}
## Get basic information
## format the log-normalized gene expression based on the selected genes.




setClass("PRECASTObj", slots=list(
  seuList = "ANY",
  seulist = "ANY",
  AdjList = "ANY", 
  parameterList= "list", 
  resList = "ANY",
  project = "character"
) )



CreatePRECASTObject <- function(seuList,  project = "PRECAST", numCores_sparkx=1, 
                                   gene.number=2000,customGenelist=NULL, premin.spots = 20,  
                                   premin.features=20, postmin.spots=15, postmin.features=15,
                              rawData.preserve=FALSE,verbose=TRUE){
  
  # premin.spots = 20;  premin.features=20; postmin.spots=15; postmin.features=15;verbose=TRUE
  #suppressMessages(require(Seurat))
  

  
  
  ## inheriting
  object <- new(
    Class = "PRECASTObj",
    seuList = seuList,
    seulist = NULL,
    AdjList = NULL, 
    parameterList= list(),
    resList = NULL,
    project = project
  )
  seuList <- object@seuList 
  if(verbose)
    message("Filter spots and features from Raw count data...")
  seuList <- lapply(seuList, filter_spot, premin.features)
  seuList <- pbapply::pblapply(seuList, filter_gene, premin.spots)
  message(" \n ")
  
  if(is.null(customGenelist)){
    seuList <- pbapply::pblapply(seuList ,DR.SC::FindSVGs, nfeatures=gene.number, 
                                 num_core=numCores_sparkx, verbose=verbose)
    spaFeatureList <- lapply(seuList, DR.SC::topSVGs, ntop = gene.number)
    spaFeatureList <- lapply(spaFeatureList, function(x) x[!is.na(x)])
    if(any(sapply(spaFeatureList, length)< gene.number)){
      gene.number_old <- gene.number
      gene.number <- min(sapply(spaFeatureList, length))
      warning(paste0("Number of SVGs in one of sample is less than ", gene.number_old, ", so set minimum number of SVGs as gene.number=", gene.number) )
    }
    if(verbose)
      message("Select common SVGs  for multiple samples...")
    genelist <- selectIntFeatures(seuList, spaFeatureList=spaFeatureList, IntFeatures=gene.number)
  }else{
    genelist <- customGenelist
  }
  
  seulist <- lapply(seuList, function(x) x[genelist, ])
  if(verbose)
     message("Filter spots and features from SVGs count data...")
  seulist <- lapply(seulist, filter_spot, postmin.features)
  seulist <- pbapply::pblapply(seulist, filter_gene, postmin.spots)
  seulist <- lapply(seulist, NormalizeData, verbose=verbose)
  object@seulist <- seulist
  
  if(!rawData.preserve){
    object@seuList <- NULL
  }
  return(object)
}


AddAdjList <- function(PRECASTObj, type="fixed_distance", platform="Visium", ...){
  if(is.null(PRECASTObj@seulist)) stop("The slot seulist in PRECASTObj is NULL!")
  
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  if(tolower(type)=='fixed_distance'){
    if(tolower(platform) %in% c("st", "visium")){
      AdjList <- pbapply::pblapply(posList, getAdj_reg, platform=platform)
    }else{
      AdjList <- pbapply::pblapply(posList, function(x, ...)getAdj_auto(x, ...))
    }
  }else if (tolower(type) == "fixed_number") {
    AdjList <- pbapply::pblapply(posList, getAdj_fixedNumber, ...)
  } else {
    stop("AddAdjList: Unsupported type \"", type, "\".")
  }
  
  # AdjList <- pbapply::pblapply(posList, function(x)getAdj_auto(x))
  PRECASTObj@AdjList <- AdjList
  return(PRECASTObj)
}

AddParSetting <- function(PRECASTObj, ...){
  PRECASTObj@parameterList <- model_set(...)
  return(PRECASTObj)
}



PRECAST <- function(PRECASTObj, K=NULL, q= 15){
  # suppressMessages(require(Matrix))
  # suppressMessages(require(Seurat))
  if(is.null(K)) K <- 4:12
  
  if(is.null(PRECASTObj@seulist)) stop("The slot seulist in PRECASTObj is NULL!")
  
  XList <- lapply(PRECASTObj@seulist, function(x) Matrix::t(x[["RNA"]]@data))
  XList <- lapply(XList, scale, scale=FALSE)
  PRECASTObj@resList <- ICM.EM_structure(XList, K=K, q=q, AdjList = PRECASTObj@AdjList, 
                             parameterList = PRECASTObj@parameterList)
  
  return(PRECASTObj)
}

## select model
selectModel.PRECASTObj <- function(obj, criteria = 'MBIC',pen_const=1, return_para_est=FALSE){
  reslist <- selectModel.SeqK_PRECAST_Object(obj@resList, pen_const = pen_const, criteria = criteria, return_para_est)
  obj@resList <- reslist
  return(obj)
}


# load("./data/Housekeeping_Genes_Mouse.RData")
# head(Mouse_HK_genes)
# load("./data/Housekeeping_GenesHuman.RData")
# Human_HK_genes <- Housekeeping_Genes
# head(Human_HK_genes)
# colnames(Human_HK_genes)[1:2]<- colnames(Mouse_HK_genes)[2:1] <- c('Ensembl', 'Gene')
# usethis::use_data(Mouse_HK_genes, Human_HK_genes, overwrite = T)


get_correct_exp <- function(XList, RfList, housekeep, q_unwanted=10){
  
  
  if(!all(sapply(XList, is.matrix))){
    XList <- lapply(XList, as.matrix)
  }
  MList <- pbapply::pblapply(XList, function(x) wpca(x[,housekeep], q=q_unwanted, F)$PCs)
  M0 <- matlist2mat(MList)
  rm(MList)
  Rf <- matlist2mat(RfList)
  rm(RfList)
  XList <-  lapply(XList, scale, scale=FALSE)
  X0 <- matlist2mat(XList)
  rm(XList)
  lm1 <- lm(X0~ 0+ cbind(Rf, M0))
  hK <- ncol(Rf)
  coefmat <- coef(lm1)[-c(1:hK),]
  rm(lm1)
  hX <- X0 - M0 %*% coefmat
  return(hX)
}

get_correct_mean_exp <- function(XList,  hVList, hW){
  
  XList <- lapply(XList, scale, scale=FALSE)
  r_max <- length(XList)
  X0 <- XList[[1]]
  hV0 <- hVList[[1]]
  if(r_max>1){
    for(r in 2:r_max){
      X0 <- rbind(X0, XList[[r]])
      hV0 <- rbind(hV0, hVList[[r]])
      
    }
  }
  
  X0 - hV0%*% base::t(hW)
  
}

IntegrateSpaData <- function(PRECASTObj, species="Human", custom_housekeep=NULL){
  # suppressMessages(require(Matrix))
  # suppressMessages(require(Seurat))
  if(is.null(PRECASTObj@seulist)) stop("The slot seulist in PRECASTObj is NULL!")
  
  XList <- lapply(PRECASTObj@seulist, function(x) Matrix::t(x[["RNA"]]@data))
  n_r <- length(XList)
  for(r in 1:n_r){
    colnames(XList[[r]]) <- firstup(colnames(XList[[r]]))
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
       intersect((genelist), Mouse_HK_genes$Gene)
      },
    mouse={
      #data(Mouse_HK_genes)
      intersect((genelist), Mouse_HK_genes$Gene)
    },
    unkown={
      character()
    }
  )
  houseKeep <- c(houseKeep, custom_housekeep)
  if(length(houseKeep) < 5){
    message("Using only PRECAST results to obtain the batch corrected gene expressions
 since species is unkown or the genelist in PRECASTObj has less than 5 overlapp 
 with the housekeeping genes of given species.")
    message("Users can specify the custom_housekeep by themselves to use the housekeeping genes based methods.")
    hX <- get_correct_mean_exp(XList,PRECASTObj@resList$hV, PRECASTObj@resList$hW)
  }else{
    message("Using bouth housekeeping gene and PRECAST results to obtain the batch corrected gene expressions.")
    hX <- get_correct_exp(XList, PRECASTObj@resList$Rf, houseKeep, q_unwanted=min(10, length(houseKeep)))
  }
  meta_data <- data.frame(batch=factor(get_sampleID(XList)), cluster= factor(unlist(PRECASTObj@resList$cluster)))
  row.names(meta_data) <- row.names(hX)
  rm(XList)
  count <- sparseMatrix(i=1,j=1, x=0, dims=dim(t(hX)))
  row.names(count) <- colnames(hX)
  colnames(count) <- row.names(hX)
  seuInt <- CreateSeuratObject(counts = count, assay = 'PRE_CAST', meta.data=meta_data)
  seuInt[['PRE_CAST']]@data <- t(hX)
  rm(hX)
  
  # seuInt <- CreateSeuratObject(assay, meta.data=meta_data, assay = 'PRECAST')
  seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$hZ), seuInt, embed_name = 'PRECAST', assay='PRE_CAST')
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  seuInt <- Add_embed(matlist2mat(posList), seuInt, embed_name = 'position', assay='PRE_CAST')
  Idents(seuInt) <- factor(meta_data$cluster)
  return(seuInt)
}


# plot cluster spatial heatmap, or tSNE RGB plot, or UMAP RGB plot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

SpaPlot <- function(seuInt, batch=NULL, item=NULL, point_size=2,text_size=16, 
                    cols=NULL,font_family='', border_col="gray10",
                    fill_col='white', ncol=2, combine = TRUE, title_name="Sample"){
  if(is.null(batch)){
    batch_vec <- unique(seuInt$batch)
  }else{
    batch_vec <- batch
  }
  if(is.null(item)) item <- "ident"
  if(is.null(cols)&& item != "RGB_UMAP" && item!="RBG_tSNE"){
    # to determine the number of colors
    meta_data <- seuInt@meta.data
    meta_data$ident <- Idents(seuInt)
    nn <- max(length(unique( meta_data$ident)), length(unique(meta_data[,item])))
    cols <- gg_color_hue(nn)
  }
  
  
  
  pList <- list()
  k <- 1
  item_null_flag <- FALSE
  for(batchi in batch_vec){
    # batchi <- (batch_vec[2])
    seu <- subset(seuInt, batch==batchi)
    meta_data <- seu@meta.data
    meta_data$ident <- Idents(seu)
     
    embed_use <- seu@reductions$position@cell.embeddings
    if(item %in% colnames(meta_data)){
      sort_id <- sort(as.numeric(as.character((unique(meta_data[, item])))))
      p1 <- plot_scatter(embed_use, meta_data, label_name=item, 
                         point_size=point_size,palette_use =cols[sort_id])
    }else if(item=="RGB_UMAP"){
      p1 <- plot_RGB(embed_use, seu@reductions$UMAP3@cell.embeddings, pointsize = point_size)
    }else if(item=="RGB_TSNE"){
      p1 <- plot_RGB(embed_use, seu@reductions$tSNE3@cell.embeddings, pointsize = point_size)
    }else{
      stop("SpaPlot: check the value of argument:item")
    }
    p1 <- p1 + mytheme_graybox(base_size = text_size, base_family = font_family, bg_fill = fill_col,
                          border_color = border_col) 
    if(!is.null(title_name)){
      p1 <- p1 + ggtitle(label=paste0(title_name, batchi))
    }
    
    pList[[k]] <- p1
    k <- k + 1
    if(item_null_flag){
      item <- NULL
    }
  }
  if(combine){
    pList <- cowplot::plot_grid(plotlist = pList, ncol=ncol)
  }
  return(pList)
}
dimPlot <- function(seuInt, item=NULL, reduction=NULL, point_size=1,text_size=16, 
                    cols=NULL,font_family='', border_col="gray10",
                    fill_col="white"){
  
  
  
  
  meta_data <- seuInt@meta.data
  if(is.null(item)){
    meta_data$ident <- Idents(seuInt)
    item <- "ident"
    
  }
  if(is.null(cols)){
    # to determine the number of colors
    nn <- length(unique( meta_data[,item]))
    cols <- gg_color_hue(nn)
  }
  if(is.null(reduction)){
    n_re <- length(seuInt@reductions)
    reduction <- names(seuInt@reductions)[n_re]
  }
  
  
  
  
  
  embed_use <- seuInt[[reduction]]@cell.embeddings[,c(1,2)]
  p1 <- plot_scatter(embed_use, meta_data, label_name=item, 
                     point_size=point_size,palette_use =cols)
 
  p1 <- p1 + mytheme_graybox(base_size = text_size, base_family = font_family, bg_fill  = fill_col,
                          border_color = border_col)
 
  return(p1)
}




# seuList <- gendata_seulist()
# PRECASTObj <-  CreatePRECASTObject(seuList)
# PRECASTObj <-  AddAdjList(PRECASTObj)
# PRECASTObj <- AddParSetting(PRECASTObj, maxIter=3)
# PRECASTObj <- PRECAST(PRECASTObj, K=4:5)
# resList <- PRECASTObj@resList
# PRECASTObj@resList <- resList
# PRECASTObj <- selectModel.PRECASTObj(PRECASTObj)

# usethis::use_data(seuInt, overwrite = T)
# seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
#seuInt

# p12 <- SpaPlot(seuInt, batch=NULL,point_size=2, combine=T)
# ggsave(file='tmp.png',plot=p12, width=9)
# seuInt <- AddUMAP(seuInt)
# SpaPlot(seuInt, batch=NULL,item='RGB_UMAP',point_size=2, combine=T, text_size=15)
# seuInt <- AddTSNE(seuInt, 2)

# dimPlot(seuInt,  font_family='serif') # Times New Roman
#  dimPlot(seuInt, reduction = 'UMAP3', item='cluster', font_family='serif') 

# head(seuInt@reductions$tSNE2@cell.embeddings)
# head(seuInt@reductions$position@cell.embeddings)
# DimPlot(seuInt, reduction = 'position')
# DimPlot(seuInt, reduction = 'tSNE')
# DimPlot(seuInt, reduction = 'PRECAST')
# library(Seurat)
# data("pbmc_small")
# AddUMAP(pbmc_small, reduction = 'pca')
# pbmc_small[["PRECAST"]] <- pbmc_small@assays$RNA
# pbmc_small <- AddUMAP(pbmc_small, reduction = 'pca', n_comp = 2)
# DefaultAssay(pbmc_small) <- "PRECAST"
# head(pbmc_small)
# dimPlot(pbmc_small, reduction = 'UMAP', item='RNA_snn_res.0.8')
# DEG analysis ------------------------------------------------------------

# dat_deg <- FindAllMarkers(seuInt)
# library(dplyr)
# n <- 10
# dat_deg %>%
#   group_by(cluster) %>%
#   top_n(n = n, wt = avg_log2FC) -> top10
# 
# seuInt <- ScaleData(seuInt)
# seus <- subset(seuInt, downsample = 400)
# color_id <- as.numeric(levels(Idents(seus)))
# 
# 
# 
# ## HeatMap
# p1 <- doHeatmap(seus, features = top10$gene, cell_label= "Domain",
#                 grp_label = F, grp_color = NULL,
#                 pt_size=6,slot = 'scale.data') + 
#   theme(legend.text = element_text(size=16),
#         legend.title = element_text(size=18, face='bold'),
#         axis.text.y = element_text(size=7, face= "italic", family='serif'))
