___
### phyloDPM 
#### Function for analysis of microbiome data from 16S rRNA gene sequencing.
##### Version: 1.0
___
*Author: Denis Awany <<awanydenis@gmail.com>>*


The R function `PhyloDPM()` implemements a Bayesian statistical method I developed for analysis of microbiome data from 16S rRNA sequencing. The framework, based on the Dirichlet Process Random (DP) Effects Model, is used to identify microbial taxa (Operational Taxonomic Units (OTUs) asscoiated with a given host phenotype. 

The method uses a weighted combination of phylogenetic- and radial basis function (RBF)-defined kernels to model microbial taxa effects, and a non-parametrically defined latent variable (model as a Dirichlet Process) to model latent heterogeneity across samples.

Given the complexity of the dependency structure of the model imposed by the DP component, which make the full conditionals for the parameters to assume awkward forms, direct Markov Chain Monte Carlo (MCMC) sampling using the Gibbs sampler is impossible. To circumvent this, posterior distribution of the parameters are here estimated by implementing what can be thought of as a `slice-Gibbs-within-Metropolis` algorithm.

The function requires, at minimum, as input: a vector of phenotypes, a matrix of microbiome abundance, and phylogenetic tree specified as an object of class matrix or phylo.


##### Implement in R by running:

```R
model <- phyloDPM(pheno=y, otu.data=X, mTree=mTree, Nsim=10000, thin=1, burnin=NULL, post.plots=FALSE, mu=0, alpha.par=NULL,
 tau2.par=c(20,20), sig2.bet.par=c(20,15), rbf.par=0.001, saveSim=FALSE, evol.rate=10^3, otu.corr=TRUE,...) 
```
where `y` is a vector of binary response, `X` is a `n` by `p` matrix/dataframe of `n` samples and `p` taxa.
The phylogenetic tree `mTree` can be created using [qiime](https://qiime2.org/), and pruned to a given taxonomic level (such as genus, specie, phylum) using tools such as [picante](https://doi.org/10.1186/s40168-023-01488-z). 

See a complete example of implementation in the script `simulation_analysis.R`.
