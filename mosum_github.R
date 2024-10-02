#' Multiple change point detection under a static factor model
#'
#' Produces estimates of multiple change points due to changes 
#' in the loadings under a static factor model.
#'
#' See Barigozzi, Cho & Trapani (2024) for further details.
#'
#' @param x data matrix with the rows containing the variables
#' @param G bandwidth
#' @param r factor number under the observationally equivalent factor model with 
#' time-invariant loadings; if \code{r = NULL}, an estimate is produced as described
#' in Alessi, Barigozzi & Capasso (2010).
#' @param V.diag boolean; if \code{V.diag = TRUE}, only the diagonal entries of
#' the estimated (long-run) covariance matrix are used for standiardisation of
#' MOSUM statistics
#' @param lrv boolean; if \code{lrv = TRUE}, the long-run covariance matrix of 
#' vech of the outer product of estimated factors are used for standiardisation of
#' MOSUM statistics; if \code{lrv = FALSE}, the variance matrix is used
#' @param m integer; bandwidth for the long-run covariance matrix estimation. 
#' The default choice is \code{floor(dim(x)[2]^.25)}
#' @param sig.lev significant level. The default choice is \code{0.05}
#' @param thr.max boolean; when \code{r} is too large, the asymptotic null 
#' distribution may return a negative critical value. To prevent this, we take 
#' the maximum of the critical values from the null distributions associated with
#' a range of \code{d = 1, ...,  r (r + 1)/2}
#' @param thr.adj determines the exponent of the multiplicative factor 
#' \code{log(dim(x)[2]/G)^a} multiplied to the critical value from the null distribution
#' @param eta a positive numeric value for the minimal mutual distance of change 
#' point estimates relative to \code{G}; the default choice is \code{eta = 0.6}
#' @return a list containing the following fields: 
#' \item{cpts}{a vector containing the estimated change point locations}
#' \item{stat}{a series of MOSUM statistic values}
#' \item{thr}{threshold}
#' \item{G}{bandwidth}
#' \item{flag}{boolean; whether the (long-run) covariance matrix has to be 
#' modified for its inversion}
#' @examples
#' \donttest{
#' set.seed(1234)
#' dd <- duan_dgp(n = 1000, p = 500, type = c(2, 1, 4), dep = TRUE)
#' G0 <- round(dim(dd$x)[2]^(max(2 / 5, 1 - min(1, 
#' log(dim(dd$x)[1])/log(dim(dd$x)[2])))) * log(dim(dd$x)[2])^1.1)
#' out <- mosum.fts(x = dd$x, G = G0)
#' out$cpts
#' plot(out$stat, type = 'l')
#' abline(h = out$thr, col = 3)
#' abline(v = out$cpts, col = 2); abline(v = dd$k0, col = 4)
#' }
#' @references Barigozzi, M., Cho, H., & Trapani, L. (2024) Moving sum procedure 
#' for multiple change point detection in large factor models.
#' @references Alessi, L., Barigozzi, M., & Capasso, M. (2010) Improved 
#' penalization for determining the number of factors in approximate static 
#' factor models. Statistics and Probability Letters, 80:1806-1813.
#' @export
mosum.fts <- function(x, G, r = NULL, 
                      V.diag = TRUE, lrv = TRUE, m = floor(dim(x)[2]^.25),
                      sig.lev = .05, thr.max = TRUE, thr.adj = .2, eta = .6){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  flag <- FALSE
  
  if(is.null(r)) r <- median(abc.factor.number(x)$r[4:6])
  d <- r * (r + 1) / 2
  
  sv <- svd(x, nu = 0, nv = r)
  g <- t(sv$v) * sqrt(n)
  ind <- lower.tri(diag(1, r), diag = TRUE)
  gg <- matrix(0, nrow = d, ncol = n)
  for(tt in 1:n) gg[, tt] <- c((g[, tt] %o% g[, tt])[ind])
  
  V <- (gg - c(diag(1, r)[ind])) %*% t(gg - c(diag(1, r)[ind])) / n
  if(m >= 1 && lrv){
    for(ll in 1:m){
      tmp <-  (gg[, 1:(n - ll), drop = FALSE] - c(diag(1, r)[ind])) %*% t(gg[, 1:(n - ll) + ll, drop = FALSE] - c(diag(1, r)[ind])) / n
      V <- V + (1 - ll/m) * (tmp + t(tmp))
    }
  }
  
  if(V.diag){
    dd <- diag(V)
    if(min(dd) < 0){
      dd <- dd * (dd >= min(dd[dd > 0])) + min(dd[dd > 0]) * (dd < min(dd[dd > 0]))
      flag <- TRUE
    }
    Vsq <- diag(1 / sqrt(dd))
  } else{
    eig <- eigen(V, symmetric = TRUE)
    dd <- eig$values
    if(min(dd) < 0){
      dd <- dd * (dd >= min(dd[dd > 0])) + min(dd[dd > 0]) * (dd < min(dd[dd > 0]))
      flag <- TRUE
    }
    Vsq <- diag(1 / sqrt(dd), d) %*% t(eig$vectors)
  }
  Vgg <- Vsq %*% gg
  
  stat <- calc_mosum(Vgg, G, boundary.ext = FALSE)
  
  thr <- get_thr(n, G, d, sig.lev, thr.max) * log(n/G)^thr.adj
  exceedings <- (stat > thr)

  localMaxima <- (c((diff.default(stat) < 0), NA) & c(NA, diff.default(stat) > 0))
  p.candidates <- which(exceedings & localMaxima)
  cpts <- mosum:::eta_criterion_help(p.candidates, stat, eta, G, G)
  for(cc in cpts) if(any(!exceedings[max(cc - 5, 1):min(cc + 5, n)])) cpts <- setdiff(cpts, cc)
  
  out <- list(cpts = cpts, stat = stat, thr = thr, G = G, flag = flag)
  return(out)
  
}

## misc

#' @keywords internal
calc_mosum <- function(y, G, boundary.ext = TRUE){
  
  n <- dim(y)[2]
  
  stat <- rep(0, n)
  for(tt in 1:(n - 1)){
    if(tt <= G - 1){
      if(boundary.ext){
        if(tt == 1){
          lft <- y[, 1]
          rgt <- apply(y[, 2:(2 * G), drop = FALSE], 1, sum)
        } else{
          lft <- lft + y[, tt]
          rgt <- rgt - y[, tt]
        }
        stat[tt] <- sum((lft / tt - rgt / (2 * G - tt))^2) * tt * (2 * G - tt) / (2 * G)
      }
    } else if(tt >= n - G + 1){
      if(boundary.ext){
        lft <- lft + y[, tt]
        rgt <- rgt - y[, tt]
        stat[tt] <- sum((lft / (2 * G + tt - n) - rgt / (n - tt))^2) * (n - tt) * (2 * G + tt - n) / (2 * G)
      }
    } else{
      if(tt == G){
        lft <- apply(y[, (tt - G + 1):tt, drop = FALSE], 1, sum)
        rgt <- apply(y[, (tt - G + 1):tt + G, drop = FALSE], 1, sum)
      } else{
        lft <- lft - y[, tt - G] + y[, tt]
        rgt <- rgt + y[, tt + G] - y[, tt] 
      }
      stat[tt] <- sum((lft - rgt)^2) / (2 * G)
    }
  }
  
  sqrt(stat)
  
}

#' @keywords internal
get_thr <- function(n, G, d, sig.lev = .05, thr.max = TRUE){
  if(thr.max){
    thr <- 0
    for(dd in 1:d){
      thr <- max(thr, (get_B(n, G, dd) - log(log(1 / sqrt(1 - sig.lev)))) / get_A(n, G))
    }
  } else thr <- (get_B(n, G, d) - log(log(1 / sqrt(1 - sig.lev)))) / get_A(n, G)
  thr
}

#' @keywords internal
get_pval <- function(z, n, G, d){
  1 - exp( -2 * exp(get_B(n, G, d) - get_A(n, G) * z))
}

#' @keywords internal
get_A <- function(n, G){
  sqrt(2 * log(n/G))
}

#' @keywords internal
get_B <- function(n, G, d){
  2 * log(n/G) + d/2 * log(log(n/G)) - log(2) - log(gamma(d / 2))
}

#' @importFrom graphics axis box legend par 
#' @keywords internal
abc.factor.number <- function(x, r.max = NULL, center = TRUE, 
                              p.seq = NULL, n.seq = NULL, do.plot = FALSE) {
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
  xx <- x - mean.x
  
  if(is.null(r.max)) r.max <- min(50, floor(sqrt(min(n - 1, p))))

  if(is.null(p.seq)) p.seq <- floor(4 * p / 5 + (1:10) * p / 50)
  if(is.null(n.seq)) n.seq <- floor(4 * n / 5 + (1:10) * n / 50)
  const.seq <- seq(.01, 3, by = 0.01)
  IC <- array(Inf, dim = c(r.max + 1, length(const.seq), 10, 6))
  
  for(kk in 1:min(length(n.seq), length(p.seq))) {
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    
    int <- sort(sample(n, nn, replace = FALSE))
    
    pen <- c((nn + pp) / (nn * pp) * log(nn * pp / (nn + pp)),
             (nn + pp) / (nn * pp) * log(min(nn, pp)),
             log(min(nn, pp)) / min(nn, pp))
    
    covx <- xx[, int] %*% t(xx[, int]) / nn
    
    sv <- svd(covx[1:pp, 1:pp], nu = 0, nv = 0)
    tmp <- rev(cumsum(rev(sv$d))) / pp
    if(pp > r.max) tmp <- tmp[1:(r.max + 1)]
    for (jj in 1:length(const.seq)) {
      for (ic.op in 1:3) {
        IC[1:length(tmp), jj, kk, ic.op] <-
          tmp + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
        IC[1:length(tmp), jj, kk, 3 * 1 + ic.op] <-
          log(tmp) + (1:length(tmp) - 1) * const.seq[jj] * pen[ic.op]
      }
    }
  }
  
  r.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(r.mat, c(1, 3), var)
  r.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      r.hat[ii] <- min(r.mat[max(which(ss == min(ss))),, ii]) - 1
    } else{
      if(sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
        r.hat[ii] <-
          r.mat[which(ss[-length(const.seq)] != 0 &
                        ss[-1] == 0)[1] + 1, dim(r.mat)[2], ii] - 1
      }else{
        r.hat[ii] <- min(r.mat[max(which(ss == 0)),, ii]) - 1
      }
    }
  }
  
  out <- list(r.hat = r.hat)
  attr(out, "data") <- list(Sc = Sc, const.seq = const.seq, r.mat = r.mat)
  
  if(do.plot){
    data <- attr(out, "data")
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    par(mfrow = c(2, 3))
    Sc <- data$Sc
    const.seq <- data$const.seq
    q.mat <- data$q.mat
    for(ii in 1:6){
      plot(const.seq, r.mat[, dim(r.mat)[2], ii] - 1, type = "b", pch = 1, col = 2, 
           bty = "n", axes = FALSE, xlab = "constant", ylab = "", main = paste("IC ", ii))
      box()
      axis(1, at = pretty(range(const.seq)))
      axis(2, at = pretty(range(r.mat[, dim(r.mat)[2], ii] - 1)),
           col = 2, col.ticks = 2, col.axis = 2)
      par(new = TRUE)
      plot(const.seq, Sc[, ii], col = 4, pch = 2, type = "b", 
           bty = "n", axes = FALSE, xlab = "", ylab = "")
      axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
      legend("topright", legend = c("r", "Sc"), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = "n")
    }    
  }
  
  return(out)

}

#' @importFrom stats toeplitz rnorm
#' @importFrom mvtnorm rmvnorm
#' @export
duan_dgp <- function(n, p, type, dep = FALSE){
 
  burnin <- n
  r <- r0 <- 3
  if(type[1] == 0){
    k0 <- n
    brks <- c(0, n)
  } else{
    k0 <- floor(n * 1:length(type) / (length(type) + 1))
    brks <- c(0, k0, n)
  }
  
  beta <- c(0, .3)[2]
  Omega <- toeplitz(beta^(1:p - 1))
  
  if(dep){
    rho <- .7
    alpha <- .3
  } else rho <- alpha <- 0
  
  f <- matrix(rnorm(r0 * (n + burnin)), nrow = r0)
  e <- t(mvtnorm::rmvnorm(n + burnin, sigma = Omega))
  for(tt in 2:dim(f)[2]){ 
    f[, tt] <- rho * f[, tt - 1] + f[, tt] 
    e[, tt] <- alpha * e[, tt - 1] + e[, tt] 
  }
  f <- f[, -(1:burnin), drop = FALSE]
  e <- e[, -(1:burnin)]
  
  lam <- array(0, dim = c(p, r0, length(k0) + 1))
  lam[,, 1] <- matrix(rnorm(r0 * p), nrow = p) / sqrt(r0)
  
  for(kk in 1:length(k0)){
    if(type[kk] == 1){
      C <- diag(rep(1, r0)); C[r0, r0] <- 0
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 2){
      C <- matrix(0, r0, r0)
      diag(C) <- c(.5, 1, 1.5)[1:r0]
      C[lower.tri(C, diag = FALSE)] <- rnorm(r0 * (r0 - 1)/2) 
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 3){
      m <- 1
      C <- matrix(c(1, 0, 0, 2, 1, 0, 3, 2, m), byrow = TRUE, nrow = r0)
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 4){
      lam[,, kk + 1] <- matrix(rnorm(r0 * p), nrow = p) / sqrt(r0)
      r <- r + r0
    }
  }
  
  chi <- matrix(0, nrow = p, ncol = n)
  for(jj in 1:(length(brks) - 1)){
    int <- (brks[jj] + 1):brks[jj + 1]
    chi[, int] <- lam[,, jj] %*% f[, int, drop = FALSE]
  }
  
  x <- chi + e
  
  return(list(x = x, r = r, k0 = k0, type = type))
  
}

