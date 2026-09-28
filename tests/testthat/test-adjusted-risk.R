test_that("native ridge fits agree with an independent BFGS oracle", {
  set.seed(291)
  for (case in 1:20) {
    n <- 40+case; p <- 4
    x <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1))
    y <- rbinom(n,1,plogis(x %*% c(.2,.4,-.6,.3)))
    if (case==1) y[] <- 1
    if (case==2) y[] <- 0
    w <- sample(0:3,n,TRUE); pen <- c(.1,rep(1,p-1)); lambda <- .5
    fn <- function(b) { z<-drop(x%*%b); sum(w*(pmax(z,0)-y*z+log1p(exp(-abs(z)))))+lambda/2*sum(pen*b*b) }
    gr <- function(b) drop(crossprod(x,w*(plogis(x%*%b)-y)))+lambda*pen*b
    ref <- optim(rep(0,p),fn,gr,method="BFGS",control=list(reltol=1e-14,maxit=2000))
    native <- poc_risk_fit(x,y,weights=w,lambda=lambda,penalty=pen)
    expect_equal(ref$convergence,0)
    expect_equal(native$objective,ref$value,tolerance=1e-9)
    expect_equal(unname(native$coefficients),ref$par,tolerance=2e-5)
    expect_lt(max(abs(gr(native$coefficients)))/sum(w),1e-9)
    expect_equal(native$predictions,drop(plogis(x%*%native$coefficients)),tolerance=1e-12)
    # Central differences are independent of the native gradient implementation.
    eps <- 1e-5
    numeric <- vapply(seq_len(p),function(j) { b<-native$coefficients;b[j]<-b[j]+eps
      bp<-fn(b);b[j]<-b[j]-2*eps;(bp-fn(b))/(2*eps)},0)
    expect_lt(max(abs(numeric-gr(native$coefficients))),1e-7)
  }
})

test_that("standardization preserves common weights and cumulative order", {
  set.seed(292)
  x <- cbind(1,matrix(rnorm(400),200,2)); g <- rep(0:1,each=100)
  y <- matrix(rbinom(600,1,.35),200,3)
  for (mode in c("annual","cumulative")) {
    fit <- poc_adjusted_risks(x,y,g,target_weights=seq_len(200),mode=mode)
    for (j in 1:2) for (h in 1:3)
      expect_equal(unname(fit$risks[j,h]),weighted.mean(fit$predictions[,h,j],seq_len(200)),tolerance=1e-12)
    if (mode=="cumulative") expect_true(all(apply(fit$risks,1,diff)>=0))
    expect_true(all(fit$risks>0 & fit$risks<1))
  }
  fit <- poc_adjusted_risks(x,y,g)
  duplicated_fit <- poc_adjusted_risks(x,y,g,target=rbind(x,x))
  expect_equal(fit$risks,duplicated_fit$risks,tolerance=1e-12)
  weighted <- poc_risk_fit(x,y[,1],weights=rep(2,200))
  doubled <- poc_risk_fit(rbind(x,x),rep(y[,1],2),target=x)
  expect_equal(weighted$coefficients,doubled$coefficients,tolerance=1e-10)
})

test_that("invalid input, unsupported risk sets and convergence failures are explicit", {
  x <- cbind(1,c(-1,0,1,-1,0,1)); y <- c(0,0,1,0,1,1); g <- c(0,0,0,1,1,1)
  expect_error(poc_risk_fit(x,rep(NA_real_,6)),"binary")
  expect_error(poc_risk_fit(x,y,lambda=0),"InvalidInput")
  expect_error(poc_risk_fit(x,y,weights=rep(0,6)),"InvalidInput")
  expect_error(poc_risk_fit(x,y,penalty=c(0,1)),"InvalidInput")
  expect_error(poc_risk_fit(x,y,max_iterations=1),"NoConvergence")
  expect_error(poc_risk_fit(x,y,max_iterations=1.1),"control")
  expect_error(poc_adjusted_risks(x,y,rep(0,6)),"both labels")
  expect_error(poc_adjusted_risks(x,y,g,target_weights=c(1,-1,1,1,1,1)),"weights")
  expect_error(poc_adjusted_risks(x,cbind(1,y),g,mode="cumulative"),"empty training risk set")
  expect_error(poc_adjusted_risks(x*2,y,g),"intercept")
  # Direct foreign-function calls must also enforce their contracts.
  expect_error(.Call(poc_risk_fit_zig_, 1L,1L,1L,1,2,1,1,1,1,100L,1e-10),"InvalidInput")
  expect_error(.Call(poc_risk_standardize_zig_,1L,1L,2,1,FALSE),"InvalidInput")
})
