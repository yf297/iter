dyn.load("iter")
library("GpGp")

#set.seed(1)

# data locations
nvec <- c(100,100)
locs1 <- as.matrix(expand.grid((1:nvec[1])/nvec[1],(1:nvec[2])/nvec[2]))
n1 <- 5000
locs1 <-locs1[sample(1:prod(nvec),n1),]

# prediction location
n2 <- 1
locs2 <- matrix( c(0.765,0.765), 1, 2 )

# order data locations by distance to prediction location
ord_near <- order(fields::rdist( locs1, locs2))
locs1 <- locs1[ord_near,]

# combine into one mat
locs <- rbind(locs1,locs2)


covparms <- c(1,0.1,0.1)
covmat <- exponential_isotropic( covparms, locs )
i1 <- 1:n1
i2 <- (n1+1):(n1+n2)

# extract A and b
A <- covmat[i1,i1,drop=FALSE]
b <- covmat[i1,i2,drop=FALSE]

n = n1
max_iter = 100
tol = 1e-6
r_norm = rep(0.0, max_iter)


system.time(a1 <- .C("cg", 
   A = as.double(A),
   b = as.double(b),
   n_p = as.integer(n),
   max_iter_p = as.integer(max_iter),
   tol_p = as.double(tol), 
   r_norm = as.double(r_norm)
))

system.time(a2 <- .C("scg", 
   A = as.double(A),
   b = as.double(b),
   n_p = as.integer(n),
   max_iter_p = as.integer(max_iter),
   tol_p = as.double(tol), 
   r_norm = as.double(r_norm)
))



pdf("pcg_subset.pdf", width=8, height=8)
plot(log10(a1$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm",col = "blue"  )
points(log10(a2$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = "red"  )
dev.off()

