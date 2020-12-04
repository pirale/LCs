library(reshape2)
library(igraph)
library(parallel)

#Load metadata
metadata = read.csv('metadata.csv', as.is=T)
rownames(metadata) = metadata$vesperID

#Load lists of same-patient clonal isolates for exclusion
eschDuplicates = scan('escherichiaColiSamePatientClonal.txt', what = 'character')
klebDuplicates = scan('klebsiellaPneumoniaeSamePatientClonal.txt', what = 'character')

exclude = unique(c(eschDuplicates, klebDuplicates))
metadata = subset(metadata, !vesperID %in% exclude)

speciesNames = c("Enterobacter_hormaechei", "Escherichia_coli", "Klebsiella_pneumoniae")

#Store cutoffs here
cutoffs = list()
#Store clustering here
clustering = list()
#Store scores here
allScores = NULL
for (species in speciesNames)
{
  #Load SNP matrix and exclude same patient isolates
  SNPMatrix  = read.table(paste(species, '.snpMatrix.tsv', sep = ''), header = T, as.is = T)
  names = SNPMatrix[,1]
  SNPMatrix = as.matrix(SNPMatrix[,-1])
  SNPMatrix = SNPMatrix[!names %in% exclude,!names %in% exclude]
  rownames(SNPMatrix) = colnames(SNPMatrix) = names[!names %in% exclude]
  
  #Create vector of unique SNP values
  distances = unique(as.numeric(SNPMatrix))
  distances = as.list(sort(distances[distances>=0]))
  
  SNPMatrix = melt(SNPMatrix,as.is = T)
  
  #Create clusterings with each possible SNP value
  tmpClusters = mclapply(distances, function(d)
  {
    #Create graph where edges link strains at most d SNPs apart
    edges = subset(SNPMatrix, value <= d)[,1:2]
    g = graph_from_edgelist(as.matrix(edges), directed=F)
    #Determine connected components of graph
    ccs = components(g)
    #Extract cluster membership identifiers
    clusters = sort(ccs$membership)
    clusters = clusters[clusters %in% which(ccs$csize>1)]
    
    #Convert cluster membership identifiers to list of assembly IDs and return
    if (length(clusters) > 0)
    {
      clusters = apply(t(unique(clusters)), 2, function(x) names(which(clusters == x)))
      if (class(clusters) != 'list')
      {
        tmp = vector(mode = 'list', ncol(clusters))
        for (c in 1:ncol(clusters))
          tmp[[c]] = clusters[,c]
        clusters = tmp
      }
      return(clusters)
    }
    else return(NULL)
  }, mc.cores = 12)
  names(tmpClusters) = distances
  
  #Score clusterings
  score = mclapply(tmpClusters, function(d) sum(unlist(lapply(d, function(x) {locations = metadata[x, 'Provider.Institution']; frequency = sort(table(locations), decreasing = T); if (length(frequency) == 1) return(frequency[1]);  return(frequency[1] - sum(frequency[-1]))}))), mc.cores=12)
  #Exract all scores
  allScores = rbind(allScores, data.frame(SNPs=as.integer(names(score)),score=unlist(score), species=species))
  
  #Determine highest scoring cutoff
  cutoff = names(which.max(score))
  
  #Store clustering and cutoffs
  clustering[[species]] = tmpClusters[[cutoff]]
  cutoffs[[species]] = as.numeric(cutoff)
}

