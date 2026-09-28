# Opt-in risk estimation. Design/scaling and target choice belong to the caller.
poc_risk_fit <- function(x, y, target = x, weights = rep(1, nrow(x)),
                         lambda = 1, penalty = c(.1, rep(1, ncol(x)-1)),
                         max_iterations = 100L, tolerance = 1e-10) {
  if (!is.matrix(x) || !is.numeric(x) || !is.matrix(target) || !is.numeric(target) ||
      ncol(x) != ncol(target) || ncol(x) < 1L || nrow(x) < 1L || nrow(target) < 1L)
    stop("x and target must be nonempty numeric design matrices with matching columns")
  if (!is.numeric(y) || length(y) != nrow(x) || anyNA(y) || any(!y %in% c(0,1)))
    stop("y must be a binary numeric vector, with no missing observations")
  if (!is.numeric(weights) || length(weights) != nrow(x) ||
      !is.numeric(penalty) || length(penalty) != ncol(x)) stop("invalid weights or penalty dimensions")
  if (length(lambda)!=1L || length(max_iterations)!=1L || !is.numeric(max_iterations) ||
      !is.finite(max_iterations) || max_iterations!=floor(max_iterations) ||
      max_iterations < 1 || max_iterations > 100000 || length(tolerance)!=1L)
    stop("invalid scalar optimization control")
  z <- .Call(poc_risk_fit_zig_, nrow(x), ncol(x), nrow(target), as.double(x),
    as.double(y), as.double(weights), as.double(target), as.double(penalty),
    as.double(lambda), as.integer(max_iterations), as.double(tolerance))
  list(coefficients = stats::setNames(z[3+seq_len(ncol(x))], colnames(x)),
       predictions = z[3+ncol(x)+seq_len(nrow(target))],
       objective = z[1], iterations = z[2], gradient = z[3],
       backend = "zig", lambda = lambda, penalty = penalty)
}

poc_adjusted_risks <- function(x, y, group, target = x,
                               weights = rep(1, nrow(x)),
                               target_weights = rep(1, nrow(target)),
                               mode = c("annual", "cumulative"), lambda = 1,
                               penalty = c(.1, rep(1, ncol(x)-1))) {
  mode <- match.arg(mode)
  if (!is.matrix(x) || !is.numeric(x) || !is.matrix(target) || !is.numeric(target) ||
      ncol(x)<1L || ncol(x)!=ncol(target) || nrow(x)<1L || nrow(target)<1L ||
      any(!is.finite(x)) || any(!is.finite(target)) ||
      any(x[,1]!=1) || any(target[,1]!=1)) stop("design matrices require finite values and a first intercept column of ones")
  if (is.numeric(y) && is.null(dim(y))) y <- matrix(y,ncol=1)
  if (!is.matrix(y) || !is.numeric(y) || nrow(y)!=nrow(x) || ncol(y)<1L ||
      anyNA(y) || any(!y %in% c(0,1))) stop("y must contain binary annual outcomes")
  if (length(group)!=nrow(x) || anyNA(group) || !setequal(unique(as.character(group)),c("0","1")))
    stop("group must contain both labels 0 and 1")
  if (!is.numeric(weights) || length(weights)!=nrow(x) || any(!is.finite(weights)) ||
      any(weights<0) || !is.numeric(target_weights) || length(target_weights)!=nrow(target) ||
      any(!is.finite(target_weights)) || any(target_weights<0) || sum(target_weights)<=0)
    stop("invalid training or target weights")
  nt <- ncol(y); nm <- nrow(target)
  risks <- matrix(NA_real_,2,nt,dimnames=list(c("0","1"),colnames(y)))
  predictions <- array(NA_real_,c(nm,nt,2),dimnames=list(NULL,colnames(y),c("0","1")))
  fits <- vector("list",2); names(fits) <- c("0","1")
  for (g in c("0","1")) {
    at_risk <- rep(TRUE,nrow(x)); pp <- matrix(NA_real_,nm,nt)
    fits[[g]] <- vector("list",nt)
    for (h in seq_len(nt)) {
      use <- as.character(group)==g & at_risk & weights>0
      if (!any(use)) stop("unsupported empty training risk set for group ",g," at horizon ",h)
      fit <- poc_risk_fit(x[use,,drop=FALSE],y[use,h],target,weights[use],lambda,penalty)
      pp[,h] <- fit$predictions
      fits[[g]][[h]] <- fit[c("coefficients","objective","iterations","gradient","backend","lambda","penalty")]
      fits[[g]][[h]]$n <- sum(use)
      if (mode=="cumulative") at_risk <- at_risk & y[,h]==0
    }
    z <- .Call(poc_risk_standardize_zig_,nm,nt,as.double(pp),as.double(target_weights),mode=="cumulative")
    risks[g,] <- z[seq_len(nt)]
    predictions[,,g] <- matrix(z[-seq_len(nt)],nm,nt)
  }
  list(risks=risks,predictions=predictions,fits=fits,mode=mode,backend="zig",
       target_n=nm,target_weight_sum=sum(target_weights),
       interpretation="Estimated standardized risks; causal identification and uncertainty require external justification")
}
