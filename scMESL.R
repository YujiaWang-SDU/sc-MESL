
#########################R packages that must be loaded##########################

#cpm
library("edgeR")
#rowVars
library("DelayedMatrixStats")
#prcomp_irlba
library("irlba")
#make.kNNG
library("loe")
#NbClust
library("NbClust")

################################gene.filtering###################################
gene.filtering <- function(data.list, original_dim, batch_size, ncores.ind, wdecay, seed)
{
  or <- list()
  cl <- parallel::makeCluster(3, outfile = "/dev/null")
  registerDoParallel(cl, cores = 3)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
    library(tensorflow)
  })
  or <- foreach(i = seq(3)) %dopar%
    {
      if (is.null(seed))
      {
        config <- list()
        config$intra_op_parallelism_threads <- ncores.ind
        config$inter_op_parallelism_threads <- ncores.ind
        session_conf <- do.call(tf$ConfigProto, config)
        sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
        keras::k_set_session(session = sess)
      } else {
        set.seed((seed+i))
        use_session_with_seed((seed+i))
      }
      
      data.tmp <- data.list[[i]]
      batch_size <-round(nrow(data.tmp)/50)
      
      x <- keras::layer_input(shape = c(original_dim))
      h <- keras::layer_dense(x, 50, kernel_constraint = keras::constraint_nonneg())
      x_decoded_mean <- keras::layer_dense(h, original_dim)
      vae <- keras::keras_model(x, x_decoded_mean)
      magrittr::`%>%`(vae,keras::compile(optimizer =   tf$contrib$opt$AdamWOptimizer(wdecay,1e-3), loss = 'mse'))
      
      
      his <- magrittr::`%>%`(vae,keras::fit(
        data.tmp, data.tmp,
        shuffle = TRUE,
        epochs = 10,
        batch_size = batch_size,
        verbose = 0
      )) 
      
      W <- keras::get_weights(keras::get_layer(vae, index = 2))[[1]]
      Wsd <- matrixStats::rowSds(W)
      Wsd[is.na(Wsd)] <- 0
      Wsd <- (Wsd-min(Wsd))/(max(Wsd)-min(Wsd))
      Wsd
    }
  
  parallel::stopCluster(cl)
  or <- matrixStats::rowMeans2(as.matrix(data.frame(or)))
  or
}

################################HVG#############################################
FS_HVG <-function(counts,k){
  library(Seurat)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = length(rownames(obj)),verbose = FALSE)
  sorted_indices <- as.numeric(obj@assays$RNA@meta.data$vf_vst_counts_rank)
  selected_indices <- order(sorted_indices)[1:k]
  selected_genes_data <- GetAssayData(object = obj)[selected_indices, ]
  selected_genes_matrix <- as.matrix(selected_genes_data)
    return(selected_genes_matrix)
}

################################HRG#############################################
FS_HRG <- function(norm,k){
  
  library(Seurat)
  library(HighlyRegionalGenes)
  
  obj=CreateSeuratObject(counts = norm)
  all.genes <- rownames(obj)
  obj= NormalizeData(obj)#log
  obj=ScaleData(obj,features = all.genes)
  obj=RunPCA(object =obj,features=all.genes)
  
  obj=FindRegionalGenes(obj,nfeatures =2000)
  
  gene_score=obj@assays$RNA@meta.data$HRG.score
  names(gene_score)=all.genes
  result=names(sort(gene_score,decreasing = TRUE))
  
  selected_genes <- result[1:k] 
  filtered_obj <- subset(obj, features = selected_genes)
  
  filtered_data <- GetAssayData(filtered_obj)
  filtered_matrix <- as.matrix(filtered_data)
  
  return(filtered_matrix)
}

################################Gini#############################################
FS_Gini <- function(expr_mat, k) {
  expr_mat <- expr_mat[Matrix::rowSums(expr_mat) > 0,]
  ginis <- apply(expr_mat, 1, reldist::gini)
  max_expr <- apply(expr_mat, 1, max)
  fit = loess(ginis~max_expr)
  outliers = abs(fit$residuals)
  outliers = outliers > quantile(outliers, probs=0.75)
  fit2 = loess(ginis[!outliers]~max_expr[!outliers])
  
  norm_ginis = rep(NA, times=length(ginis));
  norm_ginis[!outliers] = fit2$residuals;
  to_impute = which(is.na(norm_ginis))
  impute_loess <- function(i) {
    d = abs(max_expr-max_expr[i])
    closest = which(d[!outliers] == min(d[!outliers]))
    imputed = fit2$fitted[closest[1]]
    return(imputed)
  }
  fit_ginis = sapply(to_impute, impute_loess)
  norm_ginis[outliers] = ginis[outliers]-fit_ginis
  p = pnorm(norm_ginis, mean=mean(norm_ginis), sd=sd(norm_ginis), lower.tail=FALSE)
  names(p) = rownames(expr_mat)
  
  sorted_p_values <- sort(p)
  top_k_genes <- names(sorted_p_values)[1:k]
  selected_expr_mat <- expr_mat[top_k_genes, ]
  
  return(selected_expr_mat)
}
################################scDHA_FS########################################
scDHA_FS <- function(data = data, k = NULL, K = 3, n = 5000, ncores = 15L, gen_fil = T,sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-6
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > n))
  n <- ifelse(gen_fil, min(n, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 15
  
  if (gen_fil) {
    data.list <- lapply(seq(3), function(i) {
      if(!is.null(seed)) set.seed((seed+i))
      if(nrow(data) > 2000)
      {
        ind <- sample.int(nrow(data), 2000, replace = F, prob = sample.prob)
      } else {
        ind <- seq(nrow(data))
      }
      data.tmp <- as.matrix(data[ind,])
      data.tmp
    })
    
    or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)
    
    keep <- which(or > quantile(or,(1-min(n,original_dim)/original_dim)))
    
    da <- data[,keep]
    original_dim_reduce <- ncol(da)
    or.da <- da
  } else {
    da <- data
    original_dim_reduce <- ncol(da)
    or.da <- da
  }
  return(or.da)
}

################################datapreprocess###################################
FS_data <- function(data, normalize = T,log = T, dim1 = 5000,dim2 = 2000){
  
  library("scDHA")
  library("doParallel")
  library("foreach")
  library("keras")
  library("tensorflow")
  library("SingleCellExperiment")
  library("caret")
  
  #normalization
  if(normalize == T){
    data <- edgeR::cpm(data)
  }
  if(log == T){
    data <- log2(data+1)
  }
  data <- data[apply(data,1, function(x) sum(x>1) > floor(ncol(data)/20)),]
  data <- data[apply(data,1, function(x) sum(x>1) < floor(19*ncol(data)/20)),]
  data <- t(data)
  #feature selection 1
  if(dim(data)[2] > dim1)
  {
    data <- scDHA_FS(data = data,n = dim1)
    print("Feature selection 1 is performed")
  }else
  {
    print("Feature selection 1 is not performed")
  }
  #feature selection 2
  if(dim(data)[2] > dim2)
  {
    data1 <- FS_HVG(t(data),as.integer(dim2))
    print("Complete: HVG")
    data2 <- FS_HRG((2^(t(data))-1),as.integer(dim2))
    print("Complete: HRG")
    data3 <- FS_Gini(t(data),as.integer(dim2))
    print("Complete: Gini")
    
   #IE
    library("rsvd")
    library("IEntropy")
    res4 <- Get_entropy(t(data),1,0.05)
    genes4 <- res4$Gene[1:2000]
    data4<- t(data)[genes4,]
    print("Complete: IE")
    
       print("Feature selection 2 is performed")
  }else
  {
    print("Feature selection 2 is not performed")
  }
  
  data <- t(data)
  data1 <- data1
  data2 <- data2
  data3 <- data3
  data4 <- data4
  
  
  
  result <- list(data = data, 
                 data1 = data1, data2 = data2, data3 = data3,
                 data4 = data4)

  return(result)
  
}

################################IdentifyFeatures##############################

IdentifyFeatures <- function(inData = data, pFeat = 0.1){
  
  
  rm <- rowMeans(inData)
  rv <- rowVars(inData,useNames = TRUE)
  
  hind <- which(rm >= median(rm))
  my_rv <- rv[hind]
  names(my_rv) <- names(hind)
  my_vargenes <- names(head(x = sort(my_rv, decreasing = T), ceiling(pFeat*dim(inData)[1])) )
  
  topK <- names(head(x = sort(my_rv, decreasing = T), ceiling(0.05*dim(inData)[1])) )
  
  
  features <- list("VarGenes" = my_vargenes,
                   "topK" = topK
                   
  )
  return(features)
}

################################make_ensnet#########################################
make_ensnet <- function (dat_mat,  K = 20, nPCs = 10, nens=20)
{
  
  
  Ncells <- dim(dat_mat)[2]
  ngenes <- dim(dat_mat)[1]
  
  knn_ens <- Matrix(0, nrow = Ncells, ncol = Ncells, sparse = TRUE)
  sample_rng <- round(runif(n=nens, min = 0.5, max = 0.75), digit = 2)
  
  for (w in sample_rng){
    
    sampled_genes <- sample(rownames(dat_mat), size = ceiling(ngenes*w), replace = F)
    
    my_pca_log <- prcomp_irlba(t(dat_mat[sampled_genes,]), n = nPCs)
    
    
    PCs <- t(my_pca_log$x)
    cor_dist <- cor(PCs, method ='pearson')
    
    
    knn_adj <- make.kNNG(1-cor_dist, k = K, symm = F)
    
    idx <- which(knn_adj > 0)
    th <- quantile(cor_dist[idx], probs = 0.1)
    cor_dist[cor_dist>=th] <- 1
    cor_dist[cor_dist<th] <- 0
    
    
    diag(cor_dist) <- 0
    idx <- which(colSums(cor_dist) == 0)
    cor_dist[idx,] <- knn_adj[idx,]
    cor_dist[,idx] <- t(knn_adj[idx,])
    
    knn_ens <- knn_ens + cor_dist
  }
  return (knn_ens)
  
}

################################denoise#########################################
denoise <- function(data, nens = 20, knum = 30, alpha = 0.7, npc = 10){
  
  log_scdata <- data
  
  ens_net <- make_ensnet(log_scdata,  K = knum, nPCs = npc, nens=nens)
  
  knn_ensnet <- drop0(ens_net, tol = 0.5*(nens-1))
  
  diag(knn_ensnet) <- 0
  idx <- which(colSums(knn_ensnet) == 0)
  knn_ensnet[idx, ] <- ens_net[idx,]
  knn_ensnet[,idx ] <- t(ens_net[idx,])
  rm(ens_net)
  
  scale_knn <- knn_ensnet %*% Diagonal(x=1/colSums(knn_ensnet))
  scale_knn <- scale_knn %*% scale_knn
  
  eye <- Diagonal(dim(log_scdata)[2])
  #R <- (eye - (1-alpha) * scale_knn)
  R <- (eye - alpha * scale_knn)
  e <- t((1-alpha) * log_scdata)
  ssp <- solve(R,e)
  
  ccnt <- 2^ssp - 1
  ccnt <- edgeR::cpm(t(ccnt))
  
  results <- list()
  results$corrected_cpm <- ccnt
  results$eknnnet <- knn_ensnet
  
  return(results)
}
################################est_cls#######################################
est_cls <- function(data,log=T){
  
  if (log) {
    corrected_log <- log2(1+edgeR::cpm(data))
  } else {
    corrected_log <- data
  }

  myvar_genes <- IdentifyFeatures(corrected_log, pFeat = 0.1)
  var_genes <- myvar_genes$VarGenes
  my_pca_log <- prcomp_irlba(t(corrected_log[var_genes,]), n = 10)
  PCs <- t(my_pca_log$x)
  
  res <- NbClust(t(PCs), distance = "euclidean",
                 method = "ward.D", index = "rubin")
  
  return(res$Best.nc[1])
}


################################clust#######################################
clust <- function(data, best_nc, Knum = 20,log=T){
  # log_scdata <- log2(1+cpm(data))
  if (log) {
    log_scdata <- log2(1 + edgeR::cpm(data))
  } else {
    log_scdata <- data
  }
  
  Ncells <- dim(data)[2]
  myvar_genes <- IdentifyFeatures(log_scdata, pFeat=0.1)
  var_genes <- myvar_genes$VarGenes
  pca <- prcomp_irlba(t(log_scdata[var_genes,]), n = 10)
  km <- kmeans(pca$x, centers = Knum)
  
  sind <- which(table(km$cluster) == 1)
  if(length(sind) != 0){
    print("removing singleton nodes")
    clid <- setdiff(names(table(km$cluster)), names(sind))
    for(ss in names(sind)){
      single <- which(km$cluster == ss)
      mscore <- matrix(0, ncol = length(clid), nrow=1)
      ii <- 1
      for(cc in clid){
        cid <- which(km$cluster == cc)
        mscore[ii] <- cor(log_scdata[var_genes, single],
                          rowMeans2(log_scdata[var_genes, cid],useNames = TRUE))
        ii <- ii+1
      }
      km$cluster[single] <- clid[which(mscore == max(mscore))]
    }
    
  }
  
  memid <- km$cluster
  while(length(table(memid)) > best_nc){
    
    cormat <- matrix(0, nrow = Knum,  ncol = Knum)
    for (ix in 1: (Knum-1)){
      for(iy in (ix+1):Knum){
        idx <- which(memid == ix)
        idy <- which(memid == iy)
        cormat[ix, iy] <- cor(rowMeans2(log_scdata[var_genes, idx], useNames = TRUE),
                              rowMeans2(log_scdata[var_genes, idy], useNames = TRUE))
      }
    }
    
    mind <- which(max(cormat) == cormat, arr.ind = T)
    
    
    
    idx <- which(memid == mind[2])
    memid[idx] <- mind[1]
    
    Knum <- Knum - 1
    memid_new <- memid
    newid <- 1
    for(iz in names(table(memid))){
      memid_new[which(memid == iz)] <- newid
      newid <- newid + 1
    }
    memid <- memid_new
    
  }
  
  
  
  return(memid)
  
}
################################result################################

result<-function(collect){
  
  tran1 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  
  y1=collect[1,]
  x1=t(y1)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y1[j]==x1[i])
        tran1[i,j]=1
      else
        tran1[i,j]=0
    }
  }
  
  tran2 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y2=collect[2,]
  x2=t(y2)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y2[j]==x2[i])
        tran2[i,j]=1
      else
        tran2[i,j]=0
    }
  }
  
  
  tran3 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y3=collect[3,]
  x3=t(y3)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y3[j]==x3[i])
        tran3[i,j]=1
      else
        tran3[i,j]=0
    }
  }
  
  tran4 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y4=collect[4,]
  x4=t(y4)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y4[j]==x4[i])
        tran4[i,j]=1
      else
        tran4[i,j]=0
    }
  }
  
  
  tran5 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y5=collect[5,]
  x5=t(y5)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y5[j]==x5[i])
        tran5[i,j]=1
      else
        tran5[i,j]=0
    }
  }
  
  trantotal=(tran1+tran2+tran3+tran4+tran5)/5
  #View(trantotal)
  for (i in 1:ncol(collect)) {
    for (j in 1:ncol(collect)) {
      if(trantotal[i,j]>0.5){
        trantotal[i,j]=trantotal[i,j]
      }
      else{
        trantotal[i,j]=0
      }
    }
  }
  return(trantotal)
}

################################consClust################################    
consClust<-function(collect,num){
  library(SNFtool)
  
  group1=collect[1,]
  group2=collect[2,]
  group3=collect[3,]
  group4=collect[4,]
  group5=collect[5,]

  
  a1=as.numeric(max(group1))
  a2=as.numeric(max(group2))
  a3=as.numeric(max(group3))
  a4=as.numeric(max(group4))
  a5=as.numeric(max(group5))
  #Merge the results
  tt=result(collect)
  
  if(missing(num)==TRUE){
    c1=c(a1,a2,a3,a4,a5)
    c1=floor((median(c1)))
    res=spectralClustering(tt, c1)
  }
  else
    res=spectralClustering(tt, num)
  
  return(res)
}


################################scMESL################################
scMESL <- function(indata,num,npc = 10,log=F){
  
  data <- indata$data
  data1 <- indata$data1
  data2 <- indata$data2
  data3 <- indata$data3
  data4 <- indata$data4
  
  print(paste("Clustering:", 1))
  
  print("Start: estimating the number of true clusters")
  best_nc <- est_cls(data,log=F)
  print("Complete: estimating the number of true clusters")
  print("Start: initial clustering")
  cl_res <- clust(data, best_nc,log=F)
  print("Complete: initial clustering")
  
  print(paste("Clustering:", 2))
  rownames(data1) <- 1:nrow(data1)
  
  print("Start: estimating the number of true clusters")
  best_nc1 <- est_cls(data1,log=F)
  print("Complete: estimating the number of true clusters")
  print("Start: initial clustering")
  cl_res1 <- clust(data1, best_nc1,log=F)
  print("Complete: initial clustering")
  
  print(paste("Clustering:", 3))
  rownames(data2) <- 1:nrow(data2)
  
  print("Start: estimating the number of true clusters")
  best_nc2 <- est_cls(data2,log=F)
  print("Complete: estimating the number of true clusters")
  print("Start: initial clustering")
  cl_res2 <- clust(data2, best_nc2,log=F)
  print("Complete: initial clustering")
  
  print(paste("Clustering:", 4))
  rownames(data3) <- 1:nrow(data3)
  
  print("Start: estimating the number of true clusters")
  best_nc3 <- est_cls(data3,log=F)
  print("Complete: estimating the number of true clusters")
  print("Start: initial clustering")
  cl_res3 <- clust(data3, best_nc3,log=F)
  print("Complete: initial clustering")
  
  print(paste("Clustering:", 5))
  rownames(data4) <- 1:nrow(data4)
  
  print("Start: estimating the number of true clusters")
  best_nc4 <- est_cls(data4,log=F)
  print("Complete: estimating the number of true clusters")
  print("Start: initial clustering")
  cl_res4 <- clust(data4, best_nc4,log=F)
  print("Complete: initial clustering")
  
  collect <- matrix(rep(0,5*ncol(data)),5,ncol(data))
  collect[1,] <- as.numeric(cl_res)
  collect[2,] <- as.numeric(cl_res1)
  collect[3,] <- as.numeric(cl_res2)
  collect[4,] <- as.numeric(cl_res3)
  collect[5,] <- as.numeric(cl_res4)
  
  print("Start: consensus clustering")
  res <- consClust(collect,num)
  print("Complete: consensus clustering")
  
  n1=c(best_nc,best_nc1,best_nc2,best_nc3,best_nc4)
  n2=floor((median(n1)))
  output <- list()
  output$initial_clustering <-collect 
  output$clust_N1 <- n1
  output$clust_N <- n2
  output$consensus_clustering <- res
  
  
  print("Done!!")
  return(output)
}


################################evacluster(NMI、ARI)############################

evalcluster<-function(truelabel,predlabel){
  if(length(truelabel)!=length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)
  #Mutual information
  MI = 0.0
  for(idx in x_ids){
    for(idy in y_ids){
      idxOccur = which(truelabel==idx)
      idyOccur = which(predlabel==idy)
      idxyOccur = intersect(idxOccur,idyOccur)
      if(length(idxyOccur)>0){
        MI = MI + (length(idxyOccur)/total)*log2((length(idxyOccur)*total)/(length(idxOccur)*length(idyOccur)));
      }
    }
  }
  
  #Normalized Mutual information
  Hx = 0; #Entropies
  for(idx in x_ids){
    idxOccurCount = length(which(truelabel==idx));
    Hx = Hx - (idxOccurCount/total) * log2(idxOccurCount/total);
  }
  Hy = 0;#Entropies
  for(idy in y_ids){
    idyOccurCount = length(which(predlabel==idy));
    Hy = Hy - (idyOccurCount/total) * log2(idyOccurCount/total);
  }
  nmi = 2 * MI / (Hx+Hy)
  
  #(adjusted) Rand Index
  tab = table(truelabel,predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
  ari=c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2)
  
  out = c(nmi,ri,ari)
  names(out)=c("NMI","RI","ARI")
  return(out)
}





