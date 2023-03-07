dyn.load("iter")
library(GpGp)
library(pcg)
source("time.R")
library(xtable)

set.seed(1)


# data locations
nvec = c(5,5)
locs1 = as.matrix(expand.grid((1:nvec[1])/nvec[1],(1:nvec[2])/nvec[2]))
n1 = prod(nvec) 

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

scg = list()

batch = c(5, 10, 15)
steps = c(1, 2, 5)
prev  = c(1, 2, 5)

config = rev(expand.grid(batch, steps, prev))
config =  cbind(config, rep(0,27), rep(0,27))

colnames(config) = c("prev", "steps", "batch", "time", "iter")

c = 0
for(i in 1:length(batch)){
    for(j in 1:length(steps)){
	for(k in 1:length(prev)){
	    c = c + 1
	    obj <- time_scg(2, A, b, batch[i], steps[j], prev[k])
            config$time[c] = obj$time
	    config$iter[c] = obj$iter
	}
    }
}


config = config[order(config[,4]),]

print(xtable(config),include.rownames=FALSE)
time_chol(2, A, b)





#cols = c("red", "green", "pink", "orange", "yellow", "purple", "brown", "black", "gray")
#cols = rep(cols, 3)




#pdf("pcg_subset_batch5.pdf", width=8, height=8)
#plot(log10(cg$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = "blue")
#for(i in 1:9 ) points(log10(scg[[i]]$r_norm), type = "o", xlab = "iteration", ylab = "log10 residual norm", col = cols[[i]])
#legend( "bottomleft", legend = c( "steps = 1, prev = 1",
#				  "steps = 1, prev = 2",
#				  "steps = 1, prev = 5",
#	      			  "steps = 2, prev = 1",
#				  "steps = 2, prev = 2",
#				  "steps = 2, prev = 5",
#				  "steps = 5, prev = 1", 
#				  "steps = 5, prev = 2",
#				  "steps = 5, prev = 5"), pch = 1, col = cols )
#
#dev.off()


