data(kidney)
attach(kidney)

library(edgeR)
nf <- calcNormFactors(kidney.counts, method = "TMM")

library(fdrtool)
data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    norm.factors = nf)

##Save run time in repeated simulations
sort.list <- SortData(counts = kidney.counts, treatment = treatment, replic = replic,
                      sort.method = "paired", norm.factors = nf)
kidney.counts <- sort.list$counts
replic <- sort.list$replic
treatment <- sort.list$treatment
nf <- sort.list$norm.factors

probs <- CalcPvalWilcox(kidney.counts, treatment, sort.method = "paired", 
                        sorted = TRUE, norm.factors = nf, exact = FALSE)
weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    weights = weights, norm.factors = nf)

#specify exactly which genes to use in simulation
genes.select <- sample(1:nrow(kidney.counts), 500)
genes.diff <- sample(genes.select, 100)

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, genes.select = genes.select,
                    genes.diff = genes.diff, weights = weights, norm.factors = nf)




