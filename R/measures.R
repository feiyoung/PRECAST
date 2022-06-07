# ## average ilsi scores (Harmony, 2019, Nature Method)
# silho_allspots <- function(embeddings,  category){
#   
#   metadata <- data.frame(celltype=category)
#   dd <- dist(embeddings)
#   sh_scores_pro <- sapply(names(metadata), function(x) {
#     cluster::silhouette(as.numeric(as.factor(metadata[[x]])), 
#                         dd)[,3]
#   })
#   sh_scores_pro
# }
# 
# 
# ilsi_avg_scores_allspots <- function(embeddings,  category){
#   require(scPOP)
#   metadata <- data.frame(category=category)
#   lisi_scores_pro <- lisi(embeddings, meta_data = metadata, 'category')
#   lisi_scores_pro$category # return a vector
# }
# 
# F1_score_silho <- function(embeddings, celltype, sampleID){
#   require(scPOP)
#   metadata <- data.frame(celltype=celltype, sampleID = sampleID)
#   sh_scores_pro <- silhouette_width(embeddings, meta.data = metadata, c('celltype', "sampleID") )
#   sh_ct <- (1+sh_scores_pro[1])/2 # larger is better
#   sh_si <- (1+sh_scores_pro[2])/2 # smaller is better
#   f1_score <- (2* (1-sh_si)*sh_ct) / ((1-sh_si) + sh_ct)
#   return(f1_score)
# }