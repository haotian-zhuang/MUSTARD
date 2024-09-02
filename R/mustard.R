#' Trajectory-guided dimension reduction for multi-sample single-cell RNA-seq data
#'
#' @param expr The normalized and standardized gene expression matrix. Rows represent genes and columns represent cells
#' @param pseudotime The vector of user-provided pseudotime values
#' @param cellanno The vector indicating which sample each cell belongs to
#' @param interval Range of pseudotime (range of the original pseudotime by default)
#' @param r Number of components to decompose into (3 by default)
#' @param resolution Number of pseudotime values to evaluate in the temporal loading function (101 by default). It does not affect the sample or gene loading vector
#' @param smooth Smoothing parameter for RKHS norm (1e-3 by default). Larger value means smoother temporal loading functions
#' @param maxiter Maximum number of iterations (20 by default)
#' @param epsilon Convergence criteria for difference between iterations (1e-4 by default)
#'
#' @return A list including the estimated loadings and explained variances
#' \describe{
#'   \item{A.hat}{Sample loading matrix}
#'   \item{B.hat}{Gene loading matrix}
#'   \item{Phi.hat}{Temporal loading matrix}
#'   \item{time}{The pseudotime values where the temporal loading function is evaluated}
#'   \item{Lambda}{Eigenvalue vector}
#'   \item{r.square}{Variance explained by each component, which is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using each component}
#'   \item{accum.r.square}{Variance explained by the top components accumulated, which is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using the top components}
#' }
#' @export
#'
mustard <- function(expr, pseudotime, cellanno = NULL, interval = NULL, r = 3, resolution = 101, smooth = 1e-3, maxiter = 20, epsilon = 1e-4) {

  expr = expr[, names(pseudotime)]
  if(is.null(cellanno)) { cellanno = setNames(sub(':.*', '', names(pseudotime)), nm = names(pseudotime)) }
  samp = unique(cellanno)

  n = length(samp)
  p = nrow(expr)
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste0('Component', 1:r)
  colnames(A) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = samp
  rownames(B) = rownames(expr)

  # Calculate range.
  timestamps.all = sort(unique(pseudotime))
  if (is.null(interval)) interval = range(timestamps.all)

  # rescale the time to 0-1.
  pseudotime = (pseudotime - min(timestamps.all)) / (max(timestamps.all) - min(timestamps.all))
  interval = (interval - min(timestamps.all)) / (max(timestamps.all) - min(timestamps.all))

  res = NULL
  Lambda = rep(0, r)
  X = NULL
  y0 <- NULL
  Rsq <- accumRsq <- rep(0, r)

  ti <- vector(mode = 'list', length = n)
  names(ti) <- samp
  for (i in samp){
    temp = 1 + round((resolution-1) * (pseudotime[cellanno == i] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    ti[[i]] <- temp
  }

  tipos <- vector(mode = 'list', length = n)
  names(tipos) <- samp
  for (i in samp){
    keep <- ti[[i]]>0
    tipos[[i]] <- keep
    y0 <- c(y0, as.vector(t(expr[, cellanno == i][, keep])))
  }

  Lt = list()
  ind_vec <- NULL
  for (i in samp){
    Lt = c(Lt, list(pseudotime[cellanno == i]))
    ind_vec <- c(ind_vec, rep(i, sum(cellanno == i)))
  }
  names(Lt) <- samp

  tm <- unlist(unname(Lt))
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)

  for (s in 1:r){
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    print(sprintf("Calculate the %dth Component", s))

    # intialization of b
    data.unfold = NULL
    y <- NULL
    for (i in samp){
      data.unfold = cbind(data.unfold, expr[, cellanno == i])
      y <- c(y, as.vector(t(expr[, cellanno == i][, tipos[[i]]])))
    }
    # b.initials <- svd(data.unfold, nu=r, nv=r)$u
    b.initials <- irlba::irlba(data.unfold, nv=r)$u
    b.hat = b.initials[,1]
    # initialization of a
    a.hat <- rep(1,n)/sqrt(n)
    names(a.hat) <- samp
    # print(sprintf("Finish the %dth initialization", s))

    # iteratively update a,b,phi
    t <- 0
    dif <- 1
    while(t<=maxiter & dif>epsilon){
      # update phi:
      Ly = list()
      for (i in samp){
        Ly = c(Ly, list(a.hat[i]*as.numeric(b.hat %*% expr[, cellanno == i])))
      }
      names(Ly) <- samp
      gc()
      # print("start rkhs")
      phi.hat = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth = smooth)
      # print("finish rkhs")
      phi.hat = phi.hat / sqrt(sum(phi.hat^2))

      # update a:
      a.tilde <- rep(0,n)
      names(a.tilde) <- samp
      for (i in samp){
        t.temp <- tipos[[i]]
        a.tilde[i] <- b.hat %*% expr[, cellanno == i][, t.temp] %*% phi.hat[ti[[i]][t.temp]]
        a.tilde[i] <- a.tilde[i] / sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      a.new <- a.tilde / sqrt(sum(a.tilde^2))
      dif <- sum((a.hat - a.new)^2)
      a.hat <- a.new

      # update b:
      temp.num <- matrix(0,p,n)
      colnames(temp.num) <- samp
      temp.denom <- rep(0,n)
      names(temp.denom) <- samp
      for (i in samp){
        t.temp <- tipos[[i]]
        temp.num[,i] <- expr[, cellanno == i][, t.temp] %*% phi.hat[ti[[i]][t.temp]]
        temp.denom[i] <- sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      b.tilde <- as.numeric(temp.num %*% a.hat) / as.numeric(temp.denom %*% (a.hat^2))
      b.new <- b.tilde / sqrt(sum(b.tilde^2))
      dif <- max(dif, sum((b.hat - b.new)^2))
      b.hat <- b.new
      t <- t+1
      # print(t)
    }

    # calculate lambda
    x = NULL
    for (i in samp){
      t.temp = ti[[i]]
      t.temp <- t.temp[t.temp>0]
      x <- c(x,as.vector(t(a.hat[i]*b.hat %o% phi.hat[t.temp])))
    }
    X = cbind(X, x)
    l.fit = lm(y~x-1)
    lambda = as.numeric(l.fit$coefficients)
    A[,s] = a.hat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    Phi[,s] = t(phi.hat)
    Lambda[s] = lambda
    Rsq[s] <- summary(l.fit)$r.squared
    accumRsq[s] <- summary(lm(y0~X-1))$r.squared

    # update expr
    for (i in samp){
      temp <- tipos[[i]]
      expr[, cellanno == i][, which(temp)] = expr[, cellanno == i][, which(temp)] -
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s]))
    }
    print(paste0("Convergence reached at dif=", dif, ', iter=', t))
  }
  l.fit = lm(y0~X-1)
  Lambda = as.numeric(l.fit$coefficients)

  # revise the sign of Lambda
  for (r in 1:length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] <- -Lambda[r]
      A[,r] <- -A[,r]
    }
  }
  # revise the signs to make sure summation of phi is nonnegative
  sgn_phi <- sign(colSums(Phi))
  for (r in 1:ncol(Phi)){
    Phi[,r] <- sgn_phi[r]*Phi[,r]
    A[,r] <- sgn_phi[r]*A[,r]
  }
  # revise the signs to make sure summation of B is nonnegative
  sgn_B <- sign(colSums(B))
  for (r in 1:ncol(Phi)){
    B[,r] <- sgn_B[r]*B[,r]
    A[,r] <- sgn_B[r]*A[,r]
  }
  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (max(timestamps.all) - min(timestamps.all)) + min(timestamps.all)
  results = list("A.hat" = A, "B.hat" = B,
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq)
  return(results)
}


# RKHS functional regression
bernoulli_kernel <- function(x, y) {

  k1.x = x-0.5
  k1.y = y-0.5
  k2.x = 0.5*(k1.x^2-1/12)
  k2.y = 0.5*(k1.y^2-1/12)
  xy = abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy = 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy = k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

freg_rkhs <- function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth = 1e-3) {

  A <- Kmat
  for (i in names(Ly)){
    A[ind_vec == i,] <- A[ind_vec == i,]*a.hat[i]^2
  }
  cvec <- unlist(Ly)

  A.temp <- A + smooth*diag(ncol(A))

  beta <- Matrix::solve(A.temp) %*% cvec

  phi.est <- Kmat_output %*% beta
  return(phi.est)
}

#' Estimate sample loadings from the testing samples
#'
#' @param expr_test Same format as "expr" but for testing samples
#' @param pseudotime_test Same format as "pseudotime" but for testing samples
#' @param cellanno_test Same format as "cellanno" but for testing samples
#' @param res_tempted The list from \code{\link{mustard}} run on the training samples
#'
#' @return Sample loading matrix from the testing samples
#' @export
#'
est_test_sample <- function(expr_test, pseudotime_test, cellanno_test = NULL, res_tempted) {

  expr_test = expr_test[, names(pseudotime_test)]
  if(is.null(cellanno_test)) { cellanno_test = setNames(sub(':.*', '', names(pseudotime_test)), nm = names(pseudotime_test)) }
  samp_test = unique(cellanno_test)

  B <- res_tempted$B.hat
  Phi <- res_tempted$Phi.hat
  Lambda <- res_tempted$Lambda
  time_return <- res_tempted$time
  n <- length(samp_test)
  p <- nrow(B)
  r <- ncol(B)
  resolution <- length(time_return)
  A_test <- matrix(0, n, r)
  rownames(A_test) <- samp_test
  colnames(A_test) <- paste0('Component', 1:r)
  y <- NULL
  ti <- vector(mode = "list", length = n)
  names(ti) <- samp_test
  # get the coordinate of observed time points in the returned time grid
  for (i in samp_test){
    ti[[i]] <- sapply(pseudotime_test[cellanno_test == i], function(x){which.min(abs(x - time_return))})
    keep <- ti[[i]]>0
    y <- c(y, as.numeric(t(expr_test[, cellanno_test == i][, keep])))
  }
  mf.new <- expr_test

  for (s in 1:r){
    for (i in samp_test){
      t_temp <- ti[[i]]>0
      A_test[i,s] <- B[,s] %*% mf.new[, cellanno_test == i][, t_temp] %*% Phi[ti[[i]][t_temp],s]
      A_test[i,s] <- A_test[i,s] / sum((Phi[ti[[i]][t_temp],s])^2) / Lambda[s]
      mf.new[, cellanno_test == i][, t_temp] <- mf.new[, cellanno_test == i][, t_temp] -
        Lambda[s] * A_test[i,s] * (B[,s] %*% t(Phi[ti[[i]][t_temp],s]))
    }
  }
  return(A_test)
}

#' Visualize the temporal loadings
#'
#' @param res The list from \code{\link{mustard}}
#' @param r The vector of components to plot
#' @param ... Arguments in \code{ggplot2::geom_line(...)}
#'
#' @return A ggplot2 plot
#' @export
#'
plot_temporal_loading <- function(res, r = NULL, ...) {
  
  Phi.data <- res$Phi.hat
  if(is.null(r)) { r <- 1:ncol(Phi.data) }
  Phi.data <- Phi.data[, r, drop = FALSE]
  ntime <- nrow(Phi.data)
  Phi.data <- data.frame(time = res$time, value = as.vector(Phi.data),
                         component = as.factor(as.vector(t(matrix(rep(r, ntime), length(r),)))))
  ggplot(data = Phi.data, aes(x = time, y = value, color = component)) +
    geom_line(...) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Pseudotime", y = "Value", title = "Temporal loading", color = "Component")
}

#' Visualize the gene loadings
#'
#' @param res The list from \code{\link{mustard}}
#' @param i The component to plot in the x-axis
#' @param j The component to plot in the y-axis
#' @param ... Arguments in \code{ggplot2::geom_point(...)}
#'
#' @return A ggplot2 plot
#' @export
#'
plot_gene_loading <- function(res, i = 1, j = 2, ...) {
  
  B.data <- data.frame(res$B.hat)
  ggplot(data = B.data, aes(x = B.data[, i], y = B.data[, j])) +
    geom_point(...) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = paste("Component", i), y = paste("Component", j), title = "Gene loading")
}

#' Visualize the sample loadings
#'
#' @param res The list from \code{\link{mustard}}
#' @param meta The sample-level metadata
#' @param group The phenotype of interest
#' @param i The component to plot in the x-axis
#' @param j The component to plot in the y-axis
#' @param ... Arguments in \code{ggplot2::geom_point(...)}
#'
#' @return A ggplot2 plot
#' @export
#'
plot_sample_loading <- function(res, meta, group, i = 1, j = 2, ...) {
  
  A.data <- merge(data.frame(res$A.hat), meta, by = "row.names", all = F)
  rownames(A.data) <- A.data[, 1]
  A.data <- A.data[, -1]
  ggplot(data = A.data, aes(x = A.data[, i], y = A.data[, j], color = !!sym(group), shape = !!sym(group))) +
    geom_point(...) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = paste("Component", i), y = paste("Component", j), title = "Sample loading")
}

agg_median <- function(expr, pseudotime, interval) {
  
  expr_agg <- apply(interval, 1, function(i){
    cols <- (pseudotime >= i['start'] & pseudotime < i['end'])
    if(sum(cols) > 0){
      expr_fall <- expr[, cols, drop = F]
      expr_median <- apply(expr_fall, 1, median)
    } else {
      expr_median <- rep(NA, nrow(expr))
    }
  })
  rownames(expr_agg) <- rownames(expr)
  colnames(expr_agg) <- apply(interval, 1, mean)
  expr_agg <- expr_agg[, colSums(is.na(expr_agg)) == 0]
  return(expr_agg)
}

#' Cell binning
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent cells
#' @param pseudotime The vector of user-provided pseudotime values
#' @param cellanno The vector indicating which sample each cell belongs to
#' @param interval_len Number of consecutive intervals (50 by default)
#'
#' @return A list including "expr", "pseudotime", and "cellanno" of the metacells
#' @export
#'
bin_cells <- function(expr, pseudotime, cellanno, interval_len = 50) {
  
  if(is.null(cellanno)) { cellanno = setNames(sub(':.*', '', names(pseudotime)), nm = names(pseudotime)) }
  
  interval_start = seq(from = min(pseudotime), to = max(pseudotime), length.out = interval_len+1)
  interval_end = interval_start[-1]
  interval_start = interval_start[-(length(interval_start))]
  interval = cbind(interval_start, interval_end)
  colnames(interval) = c('start', 'end')
  
  datlist <- vector(mode = 'list', length = length(unique(cellanno)))
  names(datlist) <- unique(cellanno)
  for (i in unique(cellanno)) {
    expr.sample <- expr[, (cellanno == i)]
    pseudotime.sample <- pseudotime[cellanno == i]
    expr.sample.agg <- agg_median(expr = expr.sample, pseudotime = pseudotime.sample, interval = interval)
    colnames(expr.sample.agg) <- paste(i, colnames(expr.sample.agg), sep = ':')
    datlist[[i]] <- expr.sample.agg
  }
  
  expr.agg <- do.call(cbind, datlist)
  pseudotime.agg <- as.numeric(sub('.*:', '', colnames(expr.agg)))
  names(pseudotime.agg) <- colnames(expr.agg)
  cellanno.agg = setNames(sub(':.*', '', names(pseudotime.agg)), nm = names(pseudotime.agg))
  return(list(expr = expr.agg, pseudotime = pseudotime.agg, cellanno = cellanno.agg))
}

#' Gene scaling
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent cells
#' @param scale.max Maximum value to return for scaled expression (10 by default)
#'
#' @return Same format as "expr"
#' @export
#'
scale_genes <- function(expr, scale.max = 10) {
  expr_sc <- (expr - Matrix::rowMeans(expr)) / matrixStats::rowSds(expr)
  expr_sc[expr_sc > scale.max] = scale.max
  expr_sc[expr_sc < -scale.max] = -scale.max
  return(expr_sc)
}
