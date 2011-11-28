.areasum <-
function(g, n, nsize, vec){
   tmp <- .Fortran("dareachar", g=as.integer(g), n=as.integer(n), nsize=as.integer(nsize), typ=as.integer(0), vec=as.matrix(vec), res=as.matrix(rep(1.1, g)))
   tmp$res
}

.asymhubers <-
function(y, q, k = 1.5, tol = 1e-06){
   y <- y[!is.na(y)]
   n <- length(y)
   s0 <- mad(y)
   if (s0 == 0){
      return(s = 0)
   }
   th <- 2 * pnorm(k) - 1
   beta <- th + k^2 * (1 - th) - 2 * k * dnorm(k)
   for (i in 1:30) {
      yy <- pmin(pmax(- k * s0, y), k * s0)
      # make the Huber asymmetric
      at <- yy < 0
      yy[at] <- yy[at] * (1 - q) 
      at <- yy >= 0
      yy[at] <- yy[at] * q 
      yy <- yy * 2
      ss <- sum((yy)^2) / n
      s1 <- sqrt(ss/beta)
      if (abs(s0 - s1) < tol * s0){
	 break
      } 
      s0 <- s1
   }
   return(s0)
}

.center <-
function(g, n, nsize, vec){
   tmp <- .Fortran("dcenter", g=as.integer(g), n=as.integer(n), nsize=as.integer(nsize), vec=as.matrix(vec))
   tmp$vec
}

