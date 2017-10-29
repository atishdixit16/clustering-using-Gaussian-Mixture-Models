data <- read.csv('GMMdata.csv')
k <- 3
data <- data[,-1]
data <- as.matrix(data)

unknowns <- list()

normal <- function(x,mean,coVar) {
	x <- data.matrix(x)
	mean <- data.matrix(mean)
	return( as.numeric(exp(-0.5*t(x-mean)%*%solve(coVar)%*%(x-mean))) / (sqrt(2*pi*abs(det(coVar)))) )
}

#step 1: initialize
for (i in 1:k) {
	gauss <- list()
	gauss$coeff <- (1/k)
	value <- NULL
	for (j in 1:ncol(data))
		value <- c(value, mean(data[,j]))
	gauss$mu <- value
	gauss$sigma <- cov(data)
	unknowns[[i]] <- gauss
}


#step 2 : E step
gammaMatrix <- NULL
for (i in 1:nrow(data))	{
	gamma <- NULL
	for (j in 1:k) 
		gamma <- c( gamma , unknowns[[j]][[1]]*normal(data[i,],unknowns[[j]][[2]], unknowns[[j]][[3]]) )
	gamma <- gamma/sum(gamma)
	gammaMatrix <- rbind(gammaMatrix, gamma)
}

#step 3 : M step
for (i in 1:k) {
	unknowns[[i]][[1]] = ( 1/nrow(data) )*sum(gammaMatrix[,i])
	unknowns[[i]][[2]] = ( )  / ( sum(gammaMatrix[,i]) )
}


