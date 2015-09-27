impScores.hc <-
function(cghdata.regioned, dendrogram, nclusters){

######################################################################################################
# function that calculates the maximum symmetric KL-divergence (importance scores)
######################################################################################################

# if number of clusters exceeds the number of samples: error
if (nclusters > dim(cghdata.regioned$hardcalls)[2]){ stop("more clusters than samples specified!") }

# calculated number of classes used
nclass <- dim(cghdata.regioned$softcalls)[2] / dim(cghdata.regioned$hardcalls)[2]

# extract sample to cluster ids
samples.in.clusters <- cutree(dendrogram, k=nclusters)
sample.cluster.results <- cbind(1:length(samples.in.clusters), colnames(cghdata.regioned$hardcalls), samples.in.clusters)
cluster.id <- as.numeric(sample.cluster.results[,3])

# calculate mean cluster call probability profile
data.clusts <- list()
mean.clusts <- list()
for (i in 1:max(cluster.id)){
id <- c(1:length(cluster.id))[(cluster.id==i)]
data.clusts[[i]] <- cghdata.regioned$hardcalls[,id, drop=FALSE]
minCall <- min(cghdata.regioned$hardcalls)
maxCall <- max(cghdata.regioned$hardcalls)
callLevels <- seq(from=minCall, to=maxCall, by=1)
mean.clusts[[i]] <- t(apply(data.clusts[[i]], 1, function(Z, callLevels){ as.numeric(table(factor(Z, levels=callLevels))) }, callLevels=callLevels)) / dim(data.clusts[[i]])[2]
}

# calculate maximum KL-divergence between cluster pairs for all regions
max.sym.KLdiv <- numeric()
clust.dists <- NULL
for (j in 1:dim(mean.clusts[[1]])[1]){
densities <- numeric()
for (i in 1:max(cluster.id)){
densities <- cbind(densities, mean.clusts[[i]][j,])
}
precision <- as.numeric(densities)
precision <- min(precision[(precision!=0)], na.rm=TRUE)/10
        A <- KLdiv(densities, overlap = FALSE, eps = (precision)^2)
A <- A + t(A)
max.sym.KLdiv <- c(max.sym.KLdiv, max(A))
}

# generate table of importance scores
reg.imp.scores.table <- cbind(c(1:dim(cghdata.regioned$ann)[1]), cghdata.regioned$ann[,1:4], max.sym.KLdiv)
colnames(reg.imp.scores.table)[1] <- "reg.no"
reg.imp.scores.table <- data.matrix(reg.imp.scores.table)

# return table
return(reg.imp.scores.table)
}
