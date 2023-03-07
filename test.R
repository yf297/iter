dyn.load("iter")
library(GpGp)
library(pcg)
source("time.R")
library(xtable)

set.seed(1)


# data locations
nvec = c(100,100)
locs1 = as.matrix(expand.grid((1:nvec[1])/nvec[1],(1:nvec[2])/nvec[2]))
n1 = prod(nvec) 
n1
# prediction location
n2 = 1
locs2 = matrix( c(0.765,0.765), 1, 2 )

# order data locations by distance to prediction location
ord_near = order(fields::rdist( locs1, locs2))
locs1 = locs1[ord_near,]

# combine into one mat
locs = rbind(locs1,locs2)


covparms = c(1,0.1,0.1)
covmat = exponential_isotropic( covparms, locs )
i1 = 1:n1
i2 = (n1+1):(n1+n2)

# extract A and b
A = covmat[i1,i1,drop=FALSE]
b = covmat[i1,i2,drop=FALSE]



cg = time_cg(2, A, b)
cg$time
cg$iter

#scg = time_scg(2, A, b, 15, 2, 5)
#scg$time
#scg$iter

ch = time_chol(2,A,b)
ch$time
ch$iter

cols = c("blue", "red", "green")

pdf("pcg_subset.pdf", width=8, height=8)
plot(log10(cg$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = "blue", xlim = c(0, max(ch$iter, cg$iter))  )
#points(log10(scg$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = "red")
points(log10(ch$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = "green")
legend( "bottomleft", legend = c( "cg", "chol"), pch = 1, col = cols )
dev.off()


