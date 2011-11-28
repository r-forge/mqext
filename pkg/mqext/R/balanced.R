balanced <-
function(model, k, q, init=list(d=1, v=1), store=FALSE, acc=0.00001, niter=30, nonzero=TRUE){
   # preparations
   n <- model$n
   p <- model$p
   nsize <- model$nsize
   g <- model$g
   beta <- rep(1.001, p)
   iter <- 30
   x <- model$X
   y <- model$y
   v <- init$v
   d <- init$d
   # allow for more than one k
   if (is.list(k)){
      if (any(is.na(match(c("beta", "lambda", "v"), names(k))))){
	 stop("tuning constant is not correctly specified: must be\nlist with elements 'beta', 'lambda', and 'v'!")
      }
      kbeta <- k$beta
      klambda <- k$lambda
      kv <- k$v
   }else{
      kbeta <- k
      klambda <- k
      kv <- k
   }
   # allow for more than one q
   if (is.list(q)){
      if (any(is.na(match(c("beta", "lambda", "v"), names(q))))){
	 stop("quantiles are not correctly specified: must be\nlist with elements 'beta', 'lambda', and 'v'!")
      }
      qbeta <- q$beta
      qlambda <- q$lambda
      qv <- q$v
   }else{
      qbeta <- q
      qlambda <- q
      qv <- q
   }
   # center y
   y.centered <- .center(g, n, nsize, y)
   # tau
   tau <- c(beta, v, d)
   # store the iteration-specific updates
   storage <- matrix(NA, (niter + 1), (p + 2))
   # loop 
   counter <- 0
   for( i in 1:niter){
      counter <- counter + 1
      storage[i, ] <- tau
      # compute fixed effects
      tmp <- .Fortran("dbetaiter", n=as.integer(n), p=as.integer(p), g=as.integer(g), k=as.double(kbeta), xmat=as.matrix(x), yvec=as.matrix(y), v=as.double(v), d=as.double(d), nsize=as.integer(nsize), acc=as.double(acc), beta=as.matrix(tau[1:p]), iter=as.integer(iter), converged=as.integer(1), sumwgt=as.double(1.1), info=as.integer(1), dec=as.integer(0), q=as.double(qbeta))
      beta <- tmp$beta
      res <- as.numeric(y - x %*% beta)
      # compute v
      # v <- hubers(center(g, n, nsize, res), k=k, mu=0)$s^2 * (n / (g*(nsize[1]-1)))
      v <- .asymhubers(.center(g, n, nsize, res), q=qv, k=kv)^2 * (n / (g*(nsize[1]-1)))
      # compute lambda
      # lambda <- hubers(areasum(g, n, nsize, res), k=k, mu=0)$s^2 / nsize[1]^2 
      lambda <- .asymhubers(.areasum(g, n, nsize, res), q=qlambda, k=klambda)^2 / nsize[1]^2 
      # compute a 
      a <- lambda - v / nsize[1]
      # set a to zero if it is negative
      if (nonzero){
	 if (a < 0) a <- 0
      }
      # compute d (used in the subsequent iteration)
      d <- a / v
      # test whether we have to stop updating or not
      crit <- sqrt(sum((tau - c(beta, v, d))^2)/max(1e-20, sum(tau^2)))
      if (crit <= acc){
	 break 
      }else{
	 tau <- c(beta, v, d)
      }
   }
   u <- list(beta=beta, s_e=sqrt(v), s_a=sqrt(a), iter=counter)
   if (store){
      return(storage)
   }else{
      return(u)
   }
}

