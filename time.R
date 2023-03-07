time_cg = function(repetitions, A, b){
  
    n = nrow(A)
    max_iter = 500
    tol = 1e-6

    t = rep(0, repetitions)

    for(i in 1:repetitions){
        r_norm = rep(0.0, max_iter)

	t[i] = system.time(a <- .C("cg", 
	A = as.double(A),
	b = as.double(b),
	n_p = as.integer(n),
	max_iter_p = as.integer(max_iter),
	tol_p = as.double(tol), 
	r_norm = as.double(r_norm),
	k = as.integer(0)))[[3]]
    }
    
    time = mean(t[-1])
    iter = a$k
    return(list(r_norm = a$r_norm, time = time, iter = iter))

}

time_scg = function(repetitions, A, b, batch, steps, prev){
    
    n = nrow(A)
    max_iter = 500
    tol = 1e-6
    
    t = rep(0, repetitions)

    for(i in 1:repetitions){
        r_norm = rep(0.0, max_iter)

        t[i] = system.time(a <- .C("scg", 
        A = as.double(A),
        b = as.double(b),
        n_p = as.integer(n),
        max_iter_p = as.integer(max_iter),
        tol_p = as.double(tol), 
        r_norm = as.double(r_norm),
	batch = as.integer(batch),
        steps = as.integer(steps),
        prev = as.integer(prev),
	k = as.integer(0)))[[3]]
    }
    
    
    time = mean(t[-1])
    iter = a$k
    return(list(r_norm = a$r_norm, time = time, iter = iter))

}


time_chol = function(repetitions, A, b){
    
    n = nrow(A)
    max_iter = 500
    tol = 1e-6
    
    t = rep(0, repetitions)

    for(i in 1:repetitions){
        r_norm = rep(0.0, max_iter)

        t[i] = system.time(a <- .C("chol_solve", 
        A = as.double(A),
        b = as.double(b),
        n_p = as.integer(n),
        max_iter_p = as.integer(max_iter),
        tol_p = as.double(tol), 
        r_norm = as.double(r_norm),
	k = as.integer(0)))[[3]]
    }
    
    
    time = mean(t[-1])
    iter = a$k
    return(list(r_norm = a$r_norm, time = time, iter = iter))

}
