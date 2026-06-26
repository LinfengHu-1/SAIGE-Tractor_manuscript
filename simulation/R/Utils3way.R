## Utilities for 3-way admixture simulation (AFR = 0, EUR = 1, NAT = 2)
## Adapted from Tan et al Tractor-Mix simulation_Tan/Utils.R
##
## A "person" is a list with four same-length vectors:
##   hap1, hap2  -- the two haploid dosages (0/1 per site)
##   la1,  la2   -- the two per-site local-ancestry labels (0 = AFR, 1 = EUR, 2 = NAT)
## Everything else assumes that encoding.

#### gene dropping (unchanged: local ancestry travels with the haplotype) ####
GenerateChild <- function(parent1, parent2){
  genolen = length(parent1$hap1)
  idx1 = sample(c(TRUE, FALSE), size = genolen, replace = TRUE)
  childhap1 = rep(NA, genolen); childla1 = rep(NA, genolen)
  childhap1[idx1]  = parent1$hap1[idx1];  childhap1[!idx1] = parent1$hap2[!idx1]
  childla1[idx1]   = parent1$la1[idx1];   childla1[!idx1]  = parent1$la2[!idx1]

  idx2 = sample(c(TRUE, FALSE), size = genolen, replace = TRUE)
  childhap2 = rep(NA, genolen); childla2 = rep(NA, genolen)
  childhap2[idx2]  = parent2$hap1[idx2];  childhap2[!idx2] = parent2$hap2[!idx2]
  childla2[idx2]   = parent2$la1[idx2];   childla2[!idx2]  = parent2$la2[!idx2]

  list(hap1 = childhap1, hap2 = childhap2, la1 = childla1, la2 = childla2)
}


#### Balding-Nichols MAFs for 3 populations ####
## Each population has its own FST from a shared ancestral allele frequency f.
## Default FSTs (AFR=0.15, EUR=0.10, NAT=0.15) are drawn from the literature
## (Conrad et al. 2006 NatGen; 1000 Genomes 2015; Moreno-Estrada et al. 2013 PLoS Genet).
## With these values the induced pairwise FSTs are roughly
##   AFR-EUR ~ 0.15, AFR-NAT ~ 0.19, EUR-NAT ~ 0.11, in line with published pairwise FSTs
##   for continental AFR/EUR/NAT populations (Tishkoff 2009, Moreno-Estrada 2014).
GenerateMafs3 <- function(p,
                          fst_afr = 0.15, fst_eur = 0.10, fst_nat = 0.15,
                          runif1 = 0.1,   runif2 = 0.9){
  f = runif(p, runif1, runif2)
  draw <- function(fst){ rbeta(p, (1 - fst)/fst * f, (1 - fst)/fst * (1 - f)) }
  list(maf_afr = draw(fst_afr),
       maf_eur = draw(fst_eur),
       maf_nat = draw(fst_nat),
       maf_anc = f)
}


#### Independent homogeneous population ####
## `la` is a single integer ancestry label (0, 1, or 2).
GenerateHomo <- function(mafs, la){
  list(hap1 = sapply(mafs, function(maf){rbinom(1, 1, maf)}),
       hap2 = sapply(mafs, function(maf){rbinom(1, 1, maf)}),
       la1 = rep(la, length(mafs)),
       la2 = rep(la, length(mafs)))
}


#### Independent admixed 3-way individual ####
## `props` is a length-3 probability vector c(AFR, EUR, NAT).
## Each site independently draws a local ancestry on each of the two haplotypes
## from Categorical(props), then samples the allele from that population's MAF.
GenerateAdm3 <- function(mafs_afr, mafs_eur, mafs_nat, props){
  p = length(mafs_afr)
  ## la ~ Categorical(props) encoded as 0/1/2
  la1 = sample(0:2, size = p, replace = TRUE, prob = props)
  la2 = sample(0:2, size = p, replace = TRUE, prob = props)

  hap1 = rep(NA_integer_, p); hap2 = rep(NA_integer_, p)

  ## vectorised per-ancestry sampling for speed
  for (k in 0:2) {
    mafs_k = switch(as.character(k), "0" = mafs_afr, "1" = mafs_eur, "2" = mafs_nat)
    m1 = which(la1 == k)
    if (length(m1)) hap1[m1] = rbinom(length(m1), 1, mafs_k[m1])
    m2 = which(la2 == k)
    if (length(m2)) hap2[m2] = rbinom(length(m2), 1, mafs_k[m2])
  }

  list(hap1 = hap1, hap2 = hap2, la1 = la1, la2 = la2)
}


#### Writers ####
## Tractor-Mix expects a tab file with 5 meta columns + n sample columns.
WriteGenotype <- function(G, filename){
  row.names(G) = paste0("Sample", 1:nrow(G))
  write.table(
    x = cbind(CHR = 1, POS = 1:ncol(G), SNP = ".", REF = "A", ALT = "G", t(G)),
    file = filename, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE
  )
}


#### Lambda GC helper (unchanged) ####
CalcLambdaGC = function(ps, df){
  expect_stats = qchisq(ppoints(length(ps)), df = df, lower = FALSE)
  obs_stats    = qchisq((ps),                df = df, lower = FALSE)
  list(expect_stats = expect_stats, obs_stats = obs_stats)
}


#### Ancestry-specific allele counts ####
## For a given ancestry label `query` (0, 1, or 2), count alt alleles carried
## on haplotypes whose local ancestry equals `query`.
GetAncestrySpecCount <- function(ind, query){
  ((ind$la1 == query) & (ind$hap1 == 1)) * 1L +
  ((ind$la2 == query) & (ind$hap2 == 1)) * 1L
}


#### Global ancestry proportion (for population `query`) ####
GetGlobAncestry <- function(ind, query){
  mean(c(ind$la1, ind$la2) == query)
}


#### Convenience: return all three global ancestry proportions ####
GetGlobAncestry3 <- function(ind){
  la = c(ind$la1, ind$la2)
  c(AFR = mean(la == 0), EUR = mean(la == 1), NAT = mean(la == 2))
}


#### Dirichlet sampler (no extra package dependency) ####
rDirichlet <- function(n, alpha){
  k = length(alpha)
  x = matrix(rgamma(n * k, shape = alpha, rate = 1), nrow = n, byrow = TRUE)
  x / rowSums(x)
}
