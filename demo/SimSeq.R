require(SimSeq)
data(kidney)
attach(kidney)
offset <- calcNormFactors(kidney.counts, method = "TMM")

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    n1.genes = 250, n1.diff = 50, offset = offset)

##OR
sort.list <- SortData(counts = kidney.counts, treatment = treatment, replic = replic,
                      sort.method = "paired", offset = offset)
kidney.counts <- sort.list[[1]]
replic <- sort.list[[2]]
treatment <- sort.list[[3]]
offset <- sort.list[[4]]

probs <- CalcPvalWilcox(kidney.counts, treatment, sort.method = "paired", 
                        sorted = TRUE, offset = offset, exact = FALSE)
weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    n1.genes = 250, n1.diff = 50, weights = weights, offset = offset)





