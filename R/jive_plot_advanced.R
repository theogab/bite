# # get true dist
# d  <- rbind(rnorm(1000, mean=6, sd=0.01), rnorm(1000, mean=3, sd=1.5))

# # rescale them
# rs <- scales::rescale(d, to=c(10.5,11.5))
# rs1 <- re[1:1000]
# rs2 <- re[10001:2000]

# # get density coord
# de1 <- density(rs1)
# xx1 <- de1$xx
# yy1 <- de1$yy

# de2 <- density(rs2)
# xx2 <- de2$xx
# yy2 <- de2$yy


#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))


#plot.jive(phy1, my.l, regime)

foo <- function(x){
	#plot(density(x[1],x[2]), axes = FALSE, xlab=NULL, ylab=NULL, main=NULL)
	plot(density(x[1],x[2]), main="", xlab="", ylab="", xlim=c(340,360), xaxt="n", yaxt="n", axes=F)
	abline(v=350)

}

#par(mfrow=c(5,5))
#apply(dd, 1, foo)



#plot(t1, show.tip.label = FALSE)
getTreeCoords <- function(tree, xwidth=1, yheight=0.4){
	
	d = list()
	h = vcv(tree)[1] # height
	n = dim(vcv(tree))[2] # num spec

	d$xxl <- rep(h, n) + 0.5
	d$xxu <- rep(h+xwidth, n)
	
	d$yyl <- 1:n # - yheight
	d$yyu <- d$yyl + 2*yheight
	
	return(d)
}

getCoordNDC <- function(d){

	d$xxl<-grconvertX(d$xxl, from = "user", to = "ndc")
	d$xxu<-grconvertX(d$xxu, from = "user", to = "ndc")
	
	d$yyl<-grconvertX(d$yyl, from = "user", to = "ndc")
	d$yyu<-grconvertX(d$yyu, from = "user", to = "ndc")
	
	return(d)
}
#rr1 <- getTreeCoords(t1)
#rr <- getCoordNDC(rr1)


#par(fig=c(rr$xxl[i],rr$xxu[i], rr$yyl[i], rr$yyu[i]),new = TRUE, mar = rep(0, 4))
#plot(0:1, 0:1, las = 1, main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F )

# for (i in 1:10){
	# par(fig=c(rr$xxl[i],rr$xxu[i], rr$yyl[i], rr$yyu[i]), new=TRUE, mar = rep(0, 4))
	# #plot(0:1, 0:1, las = 1, main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F )
	# hist(rnorm(1000,0.5,0.1), main="", xlab="", ylab="", xaxt="n", yaxt="n", axes=F, xlim=c(0,1))
	
# }

# require(ape)
# require(TeachingDemos)

# layout( matrix(c(1,4,4,2,4,4,3,4,4), 3, 3, byrow = TRUE))

# hist(rnorm(1000,25,3),xlab='',ylab='',main='')
# hist(rnorm(1000,25,3),xlab='',ylab='',main='')
# hist(rnorm(1000,25,3),xlab='',ylab='',main='')
# #par(mar=c(5.1,4.1,4.1,2.1))
# plot(t1, show.tip.label = TRUE, label.offset = 2, x.lim=c(0, 12.5), cex=0.9)
# #tiplabels(t1$tip.label, frame = "none")

# getTreeCoords <- function(tree, offs = 1.2, xwidth=1, yheight=0.4){
	
	# d = list()
	# h = vcv(tree)[1] # height
	# n = dim(vcv(tree))[2] # num spec

	# d$xxl <- rep(h, n) + offs
	# d$xxu <- rep(h+xwidth, n)
	
	# d$yyl <- 1:n # - yheight
	# d$yyu <- d$yyl + 2*yheight
	
	# return(d)
# }
# rr1 <- getTreeCoords(t1)

# for (i in 1:10){
	# if (i < 4){
		# par(mgp=c(0, 0, 0))
		# d <- density(rnorm(100, mean=1.5, sd=0.1))
		# subplot( plot( d ,xlab='',ylab='',main='', yaxt="n",xaxt="n", bty="n", xlim=c(-2,2), col="blue"), rr1$xxl[i], rr1$yyl[i], size=c(0.5,0.5))

	# }
	# else{
		# subplot( plot( density(rnorm(100, mean=0, sd=0.5)),xlab='',ylab='',main='', yaxt="n", xaxt="n", bty="n", xlim=c(-2,2)), rr1$xxl[i], rr1$yyl[i], size=c(0.5,0.5))
	# }
# }



# subplot(polygon(d, col="red", border="gray",xlab='',ylab='',main='', xaxt="n", yaxt="n"), 11.2, 1.8, size=c(0.5,0.5))
# subplot(plot( density(rnorm(1000)),xlab='',ylab='',main='', xaxt="n", yaxt="n", axes=F), 10.5, 4, size=c(0.5,0.5))

