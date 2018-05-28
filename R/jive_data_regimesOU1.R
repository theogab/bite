#' @title Jive test dataset: selective regimes
#' @description A vector of selective regimes generated under OU1 models
#' 
#' @details This dataset includes a set of simulated trees and trait values. The parameters used for 
#' simulation of these datasets are:
#' \itemize{
#' 	\item BM. For model.mean - BM with sig.sq = 0.5, anc.state = 350; model.var - BM with sig.sq = , anc.state = 
#' 	\item OU1. For model.mean - BM with sig.sq = 0.5, anc.state = 350; model.var - OU1 with alpha = , sig.sq = , anc.state = , theta1 = 
#' 	\item BM. For model.mean - BM with sig.sq = 0.5, anc.state = 350; model.var - OUM with alpha = , sig.sq = , anc.state = , theta1 = , theta2 = 
#' }
#' 
#' @name regimesOU1
#' @usage regimesOU1
#' @format Three phylogenetic trees of 50 species in phylo format and three data matrices with 50 rows
#' @docType data
#' @author Anna Kostikova 
#' @keywords datasets, data
NULL