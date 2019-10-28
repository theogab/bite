#' @title Format traits matrix to parse into \code{\link{make_jive}} function
#' @description This function takes a dataframe with indivdual observations ans species names and returns a traits matrix to to parse into \code{\link{make_jive}} function
#' @param obs a dataframe with species names in the first column and individual trait values in the second columns
#' @export
#' @author Theo Gaboriau
#' 


format_jive_traits <- function(obs){
  traits <- lapply(unique(obs[,1]), function(x){
    obs[obs[,1] == x,2]
  })
  nb <- max(sapply(traits, length))
  traits <- t(sapply(traits, function(x){
    c(x, rep(NA, nb - length(x)))
  }))
  rownames(traits) <- unique(obs[,1])
  return(traits)
}