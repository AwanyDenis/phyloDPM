# ------------------------------------------------------------------------------
# Simulating Microbiome Abundance Data
# ------------------------------------------------------------------------------
library(dirmult) # for simulating from DM distribution.
library(GUniFrac) # contains throat data (real data) to be used to obtain realistic sim params.
library(picante) # to prune phylogenetic tree.

# Real data for simulation parameters
data(throat.otu.tab)
data(throat.tree)

n <- 1000 # number of samples to simulate
m <- 200 # number of OUs  to simulate

# Only 200 most abundant taxa were used for parameter estimation.
otu.pruned <- throat.otu.tab[order(colSums(throat.otu.tab))[1:m]]
otu.pruned <- otu.pruned[rowSums(otu.pruned) !=0,]
DM.params <- dirmult(data=otu.pruned,epsilon=10^(-4), trace=TRUE) #obtain 'pi' and 'theta'.

# simulate OTU counts
otu.counts <- simPop(J=n, # number of samples
                     n=m,  # number of OTUs
                     pi=DM.params$pi,
                     theta=DM.params$theta)$data

# Read counts sampled from negative binom with mean 5000 and dispersion 25:
#rnbinom(n=1,mu=5000,size=25); and normalize the counts.
otu.abun <- as.data.frame(apply(otu.counts,2,
                                function(x){return(x/rnbinom(n=1,mu=5000,size=25) )}) );
colnames(otu.abun)=paste0("OTU",1:ncol(otu.counts))

# ------------------------------------------------------------------------------
write.table(as.data.frame(otu.abun), file="simulatedOTUs.txt", row.names=FALSE)
# ------------------------------------------------------------------------------
