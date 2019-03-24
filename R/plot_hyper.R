plot_hyper <- function(hpf, col = "#2e86ab", border = "#000000", ...){
  
  if(hpf(0)[[2]][[1]] == "Uniform"){
    minmax <- range(runif(1e6, min = hpf(0)[[2]][[2]][1], max = hpf(0)[[2]][[2]][2]))
    main <- sprintf("Uniform with min = %s and max = %s",round(hpf(0)[[2]][[2]][1],2), round(hpf(0)[[2]][[2]][2],2))
  }
  
  if(hpf(0)[[2]][[1]] == "Gamma"){
    minmax <- range(rgamma(1e6, shape = hpf(0)[[2]][[2]][1], scale = hpf(0)[[2]][[2]][2]))
    main <- sprintf("Gamma with shape = %s and scale = %s",round(hpf(0)[[2]][[2]][1],2), round(hpf(0)[[2]][[2]][2],2))
  }
  
  if(hpf(0)[[2]][[1]] == "Normal"){
    minmax <- range(rnorm(1e6, min = hpf(0)[[2]][[2]][1], hpf(0)[[2]][[2]][2]))
    main <- sprintf("Normal with mean = %s and variance = %s",round(hpf(0)[[2]][[2]][1],2), round(hpf(0)[[2]][[2]][2],2))
  }
  
  x <- seq(minmax[1] - diff(minmax)*0.05, minmax[2] + diff(minmax)*0.05, length.out = 1e5)
  Density <- exp(sapply(x, function(s) hpf(s)[[1]]))
  
  plot(Density~x, type = "n", main = main,...)
  
  if(is.null(border)){
    polygon(x, Density, col = col, border = NA)
  } else if(is.null(col)){
    lines(x, Density, col = border)
  } else {
    polygon(x, Density, col = col, border = NA)
    lines(x, Density, col = border)
  }
}