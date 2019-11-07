plot_hyper <- function(hpf, col = "#2e86ab", border = "#000000", ...){
  
  if(hpf(0)[[2]][[1]] == "Uniform"){
    minmax <- range(runif(1e6, min = hpf(0)[[2]][[2]][1], max = hpf(0)[[2]][[2]][2]))
    main <- substitute(U(a == min, b == max),list(min = round(hpf(0)[[2]][[2]][1],2), max = round(hpf(0)[[2]][[2]][2],2)))
  }
  
  if(hpf(0)[[2]][[1]] == "Gamma"){
    minmax <- range(rgamma(1e6, shape = hpf(0)[[2]][[2]][1], scale = hpf(0)[[2]][[2]][2]))
    main <- substitute(Gamma(alpha == a, beta == b), list(a = round(hpf(0)[[2]][[2]][1],2), b = round(hpf(0)[[2]][[2]][2],2)))
  }
  
  if(hpf(0)[[2]][[1]] == "Normal"){
    minmax <- range(rnorm(1e6, mean = hpf(0)[[2]][[2]][1], sd = hpf(0)[[2]][[2]][2]))
    main <- substitute(N(mu == m, sigma^2 == v), list(m = round(hpf(0)[[2]][[2]][1],2), v = round(hpf(0)[[2]][[2]][2],2)))
  }
  
  if(hpf(0)[[2]][[1]] == "Lognormal"){
    minmax <- range(rlnorm(1e6, mean = hpf(0)[[2]][[2]][1], sd = hpf(0)[[2]][[2]][2]))
    main <- substitute(Lognormal(mu == m, sigma^2 == v), list(m = round(hpf(0)[[2]][[2]][1],2), v = round(hpf(0)[[2]][[2]][2],2)))
  }
  
  x <- seq(minmax[1] - diff(minmax)*0.2, minmax[2] + diff(minmax)*0.2, length.out = 1e5)
  Density <- exp(sapply(x, function(s) hpf(s)[[1]]))
  
  plot(Density~x, type = "n", main = main,
       xlim = c(min(x) - diff(range(x))*.1, max(x) + diff(range(x))*.1),
       ylim = c(0, max(Density) + diff(range(Density))*.1),
       ...)
  
  if(is.null(border)){
    polygon(x, Density, col = col, border = NA)
  } else if(is.null(col)){
    lines(x, Density, col = border)
  } else {
    polygon(x, Density, col = col, border = NA)
    lines(x, Density, col = border)
  }
}