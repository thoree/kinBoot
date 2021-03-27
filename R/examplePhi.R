#'
#' Simulation experiment
#'
#' We consider parametric and nonparametric bootstrap and the kinship coefficient for some pedigrees.
#'
#' @param peds list of ped objects with allele frequencies.
#' @param idlist list of ids of of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param CItype Logical
#' @param conf.level Double
#' @param seed Integer
#' @param verbose Logical
#' @param db database
#' @param case Logical. Describes experiment.
#' #'
#' @return Returns a data frame summarising the simulation.
#' @details See [kinBoot::bootPhi()]
#'
#' @export
#'
#' @examples
#' library(forrel)
#' library(ribd)
#' library(coxed) # for bca confidence intervals
#'
#' peds = list(quadHalfFirstCousins(), doubleFirstCousins(), nuclearPed(2),
#'             halfSibPed(), cousinPed(1))
#' names(peds) = c("QHFC",  "DFC", "S", "H", "FC")
#' idlist = lapply(peds, leaves)
#' # phi = unlist(lapply(peds, function(x) kinship(x, leaves(x))))

#' N = 1; B = 4 # Increase to N = 100 and B = 1000, at least
#' seed = 1729

#' # Example 1a Many SNP-s. Intended to meet assumptions well
#' n = 1000 # no of markers
#' p = rep(0.5, n)
#' freq = list()
#' for (i in 1:n)
#' freq[[i]] =  list(afreq = c("1" = p[i], "2" = 1- p[i]))
#' db = freq
#' res1 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db, case = "1000snps")
#'
#' # Example 1b NorwegianFrequency. Medium number of markers
#' db = NorwegianFrequencies
#' res2 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db, case ="Norwegian")
#'
#' # Example 1c Few markers, 8 CODIS markers
#' codis8 = c("CSF1PO", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51")
#' db = NorwegianFrequencies[codis8]
#' res3 = examplePhi(peds, idlist, N = N, B = B, seed = seed, db = db, case = "codis8")

#' # Example 1d
#' db1 = NorwegianFrequencies[c("SE33")]
#' res4 = examplePhi(peds, idlist, N = N, B = B, seed = seed,
#'                   db = db1, case = "db1")
#'
#' # Example 1e
#' n = 10 # no of markers
#' p = 0.5
#' freq = list()
#' for (i in 1:(n-1))
#'      freq[[i]] =  list(afreq = c("1" = p, "2" = 1- p))
#' p = rep(1/100,100)
#' names(p) = 1:100
#' freq[[10]] = list(afreq = p)
#' res5 = examplePhi(peds, idlist, N = N, B = B, seed = seed,
#'                 db = freq, case = "nonid")
#'
#'  res = rbind(res1, res2, res3, res4, res5)
#'  save(res, file = "res27mar.Rdata")



examplePhi <- function(peds, idlist, N = 2, B = 10,  CItype = "perc",
                       conf.level = 0.95, seed = NULL, verbose  = TRUE, db = NULL,
                       case = "foo"){
  if(verbose) cat(date())
  if(!is.null(seed))
    set.seed(seed)
  n = length(peds)
  res = matrix(ncol = 6, nrow = 2*n)
  res = res0 = data.frame()
  for(i in 1:n){
    if(verbose) cat("\n", i, " ")
    peds[[i]] = setMarkers(peds[[i]], locusAttributes = db)
    phi = ribd::kinship(peds[[i]], idlist[[i]])
    res0 =  bootPhi(ped = peds[[i]], ids = idlist[[i]], N = N, B = B,
            CItype = CItype, conf.level = conf.level)$average
    res0 = cbind(phi = rep(phi,2), res0)
    res = rbind(res, res0)
  }

  res = data.frame(case = case, par = rep(c("Y","N"), n), ped = rep(names(peds), each =2), res)
  rownames(res) = NULL
  if(verbose) cat(date(), "\n")
  res
}





