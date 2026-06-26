## 3-way admixture pedigree (adapted from simulation_Tan/Pedigree3.R)
## Same 20-person, 10-retained-descendant pedigree structure as Tan et al.
## Each founder draws a 3-way Dirichlet admixture proportion.
##
## Outputs per family:
##   GenoMatTot  -- total alt dosage (n=10 x p)
##   GenoMatAfr, GenoMatEur, GenoMatNat -- ancestry-specific alt dosages
##   LAMat_Afr, LAMat_Eur, LAMat_Nat    -- per-site local ancestry counts (0/1/2)
##   Admprop  (10 x 3)  -- per-individual global ancestry (AFR, EUR, NAT)

MakePedigree3 <- function(MAFs, dirichlet_alpha = c(1.0, 5.5, 3.5)){
  DrawFounder <- function(){
    props = as.numeric(rDirichlet(1, dirichlet_alpha))
    GenerateAdm3(MAFs$maf_afr, MAFs$maf_eur, MAFs$maf_nat, props)
  }

  Sample1  = DrawFounder()
  Sample2  = DrawFounder()
  Sample3  = DrawFounder()
  Sample4  = GenerateChild(Sample1, Sample2)
  Sample5  = GenerateChild(Sample1, Sample2)
  Sample6  = DrawFounder()
  Sample7  = GenerateChild(Sample1, Sample2)
  Sample8  = DrawFounder()
  Sample9  = GenerateChild(Sample3, Sample4)
  Sample10 = GenerateChild(Sample3, Sample4)
  Sample11 = GenerateChild(Sample3, Sample4)
  Sample12 = GenerateChild(Sample5, Sample6)
  Sample13 = GenerateChild(Sample5, Sample6)
  Sample14 = GenerateChild(Sample7, Sample8)
  Sample15 = GenerateChild(Sample7, Sample8)
  Sample16 = GenerateChild(Sample7, Sample8)
  Sample17 = DrawFounder()
  Sample18 = GenerateChild(Sample16, Sample17)
  Sample19 = GenerateChild(Sample16, Sample17)
  Sample20 = GenerateChild(Sample16, Sample17)

  kept = list(Sample9, Sample10, Sample11, Sample12, Sample13,
              Sample14, Sample15, Sample18, Sample19, Sample20)

  totRow <- function(ind){ ind$hap1 + ind$hap2 }
  laRow  <- function(ind, q){ (ind$la1 == q) * 1L + (ind$la2 == q) * 1L }

  GenoMatTot = do.call(rbind, lapply(kept, totRow))
  GenoMatAfr = do.call(rbind, lapply(kept, GetAncestrySpecCount, query = 0))
  GenoMatEur = do.call(rbind, lapply(kept, GetAncestrySpecCount, query = 1))
  GenoMatNat = do.call(rbind, lapply(kept, GetAncestrySpecCount, query = 2))
  LAMat_Afr  = do.call(rbind, lapply(kept, laRow, q = 0))
  LAMat_Eur  = do.call(rbind, lapply(kept, laRow, q = 1))
  LAMat_Nat  = do.call(rbind, lapply(kept, laRow, q = 2))

  Admprop = t(sapply(kept, GetGlobAncestry3))  ## 10 x 3

  list(GenoMatTot = GenoMatTot,
       GenoMatAfr = GenoMatAfr, GenoMatEur = GenoMatEur, GenoMatNat = GenoMatNat,
       LAMat_Afr  = LAMat_Afr,  LAMat_Eur  = LAMat_Eur,  LAMat_Nat  = LAMat_Nat,
       Admprop    = Admprop)
}


######################### MakeGenotype: both Ind and Rel ###############################
## Wrapper kept for the Power simulation step (small-n), for full-scale simulation we
## use the block-wise Pedigree3way_Ind.R / _Rel.R scripts.
MakeGenotype3 <- function(nInd, nRel, MAFs, dirichlet_alpha = c(1.0, 5.5, 3.5),
                          write = FALSE, filename){
  famsize = 10
  p = length(MAFs$maf_afr)
  if (nRel %% famsize != 0) stop("nRel must be a multiple of 10")

  ## Independents
  GTot_Ind = matrix(NA_integer_, nInd, p); GAfr_Ind = GTot_Ind
  GEur_Ind = GTot_Ind; GNat_Ind = GTot_Ind
  LAafr_Ind = GTot_Ind; LAeur_Ind = GTot_Ind; LAnat_Ind = GTot_Ind
  Admprop_Ind = matrix(NA_real_, nInd, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

  for (i in seq_len(nInd)) {
    props = as.numeric(rDirichlet(1, dirichlet_alpha))
    ind = GenerateAdm3(MAFs$maf_afr, MAFs$maf_eur, MAFs$maf_nat, props)
    GTot_Ind[i, ] = ind$hap1 + ind$hap2
    GAfr_Ind[i, ] = GetAncestrySpecCount(ind, 0)
    GEur_Ind[i, ] = GetAncestrySpecCount(ind, 1)
    GNat_Ind[i, ] = GetAncestrySpecCount(ind, 2)
    LAafr_Ind[i, ] = (ind$la1 == 0) + (ind$la2 == 0)
    LAeur_Ind[i, ] = (ind$la1 == 1) + (ind$la2 == 1)
    LAnat_Ind[i, ] = (ind$la1 == 2) + (ind$la2 == 2)
    Admprop_Ind[i, ] = GetGlobAncestry3(ind)
  }

  ## Related
  GTot_Rel = matrix(NA_integer_, nRel, p); GAfr_Rel = GTot_Rel
  GEur_Rel = GTot_Rel; GNat_Rel = GTot_Rel
  LAafr_Rel = GTot_Rel; LAeur_Rel = GTot_Rel; LAnat_Rel = GTot_Rel
  Admprop_Rel = matrix(NA_real_, nRel, 3, dimnames = list(NULL, c("AFR","EUR","NAT")))

  for (i in seq_len(nRel / famsize)) {
    fam = MakePedigree3(MAFs, dirichlet_alpha)
    rng = (famsize * (i - 1) + 1):(famsize * i)
    GTot_Rel[rng, ] = fam$GenoMatTot
    GAfr_Rel[rng, ] = fam$GenoMatAfr
    GEur_Rel[rng, ] = fam$GenoMatEur
    GNat_Rel[rng, ] = fam$GenoMatNat
    LAafr_Rel[rng, ] = fam$LAMat_Afr
    LAeur_Rel[rng, ] = fam$LAMat_Eur
    LAnat_Rel[rng, ] = fam$LAMat_Nat
    Admprop_Rel[rng, ] = fam$Admprop
  }

  GTot = rbind(GTot_Ind, GTot_Rel)
  GAfr = rbind(GAfr_Ind, GAfr_Rel)
  GEur = rbind(GEur_Ind, GEur_Rel)
  GNat = rbind(GNat_Ind, GNat_Rel)
  LAafr = rbind(LAafr_Ind, LAafr_Rel)
  LAeur = rbind(LAeur_Ind, LAeur_Rel)
  LAnat = rbind(LAnat_Ind, LAnat_Rel)
  Admprop = rbind(Admprop_Ind, Admprop_Rel)

  list(GTot = GTot, GAfr = GAfr, GEur = GEur, GNat = GNat,
       LAafr = LAafr, LAeur = LAeur, LAnat = LAnat,
       Admprop = Admprop)
}


######################### True kinship matrix from the pedigree ###############################
MakeGRM3 <- function(nInd, nRel){
  famsize = 10
  if (nRel %% famsize != 0) stop("nRel must be a multiple of 10")
  n = nInd + nRel
  GRM = matrix(0, n, n); diag(GRM) = 1

  kindf <- data.frame(
    id  = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    dad = c(0,0,0,1,1,0,1,0,3, 3, 3, 5, 5, 7, 7, 7, 0, 16,16,16),
    mom = c(0,0,0,2,2,0,2,0,4, 4, 4, 6, 6, 8, 8, 8, 0, 17,17,17),
    sex = c(1,0,1,0,1,0,1,0,0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0)
  )
  kinfam = kinship2::kinship(with(kindf, kinship2::pedigree(id, mom, dad, sex)))
  kinfam = kinfam[c(9:15, 18:20), c(9:15, 18:20)] * 2

  GRM_Rel = matrix(0, nRel, nRel)
  for (i in seq_len(nRel / famsize)) {
    rng = (famsize * (i - 1) + 1):(famsize * i)
    GRM_Rel[rng, rng] = kinfam
  }

  GRM[(nInd + 1):n, (nInd + 1):n] = GRM_Rel
  row.names(GRM) = paste0("Sample", 1:n)
  colnames(GRM) = paste0("Sample", 1:n)
  GRM
}
