sample.cluster.table <-
function(cghdata.regioned, dendrogram, nclusters){

######################################################################################################
#
# generate table of which samples belong to which cluster
#
######################################################################################################

# if number of clusters exceeds the number of samples: error
if (nclusters > dim(cghdata.regioned$hardcalls)[2]){ stop("more clusters than samples specified!") }

# generate table specifying which samples belong to which cluster
samples.in.clusters <- cutree(dendrogram, k=nclusters)
sample.cluster.results <- cbind(1:length(samples.in.clusters), colnames(cghdata.regioned$hardcalls), samples.in.clusters)
colnames(sample.cluster.results) <- c("sample.number", "sample.label", "cluster")

# return table
return(sample.cluster.results)
}
