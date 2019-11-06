beast_data <- function (file, tree) 
{
  digits <- NULL
  pp <- prop.part(tree)
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  X <- X[grep("tree TREE1[[:space:]]+=", X)]
  X <- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
  sp <- gsub("\\[(.*?)\\]", "", X)
  sp <- strsplit(sp, ",")[[1]]
  sp <- gsub("([^\\(]*)\\(", "", sp)
  sp <- unlist(sapply(strsplit(sp, ":"), function(x){
    x <- x[-2]
    x[-1] <- NA
    as.numeric(x)
  }))
  mrcas <- lapply(sp,c)
  rm <- 0
  for(i in which(is.na(sp))){
    j <- i-rm
    desc <- c(mrcas[[j-2]], mrcas[[j-1]])
    sp[i] <- which(sapply(pp, function(x) all(x%in%desc) & all(desc%in%x))) + length(tree$tip.label)
    mrcas[[j]] <- desc
    mrcas <- mrcas[-(j-1:2)]
    rm <- rm + 2
  }
  tab <- unlist(strsplit(X, "\\["))[-1]
  tab <- gsub("&|;|\\]", "", tab)
  tab <- gsub(":.+$", "", tab)
  foo <- function(x) {
    x <- unlist(strsplit(x, ","))
    x
  }
  tab <- lapply(tab, foo)
  for (i in seq(along = tab)) {
    ind <- grep("[{]", tab[[i]])
    names <- gsub("=.+$", "", tab[[i]][ind])
    tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
    tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
    tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
    tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), 
                               tab[[i]][ind + 1])
  }
  ttab <- data.frame()
  stats <- unique(gsub("=.+$", "", unlist(tab)))
  for (i in seq(along = tab)) {
    for (j in seq(along = stats)) {
      ind <- grep(paste("^", stats[j], "=", sep = ""), 
                  tab[[i]])
      if (length(ind) > 0) {
        v <- as.numeric(gsub(paste(stats[j], "=", sep = ""), 
                             "", tab[[i]][ind]))
        ttab[i, j] <- v
      }
    }
  }
  colnames(ttab) <- stats
  ttab <- cbind(nodes = sp, ttab)
}