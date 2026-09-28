# Synthetic end-to-end example; no files, credentials or acquisition.
run_pnsbib_adjusted_examples <- function() {
  library(pnsbib)
  group <- rep(0:1, each=8)
  x <- cbind(intercept=1, baseline=rep(c(-1,0,1,0,-1,0,1,0),2))
  y <- cbind(t1=rep(c(0,0,1,1,0,1,0,1),2),
             t2=rep(c(0,1,0,1,1,0,1,0),2))
  annual <- poc_adjusted_risks(x,y,group,mode="annual")
  cumulative <- poc_adjusted_risks(x,y,group,mode="cumulative")
  observed <- unclass(table(factor(group,0:1),factor(y[,1],0:1)))/length(group)
  dimnames(observed) <- list(c("0","1"),c("0","1"))
  risk <- annual$risks[,1]
  intervention <- cbind("0"=1-risk,"1"=risk)
  q <- poc_query(c("0"="0","1"="1"))
  model <- poc_model(observed,intervention,lp_backend="zig")
  complete <- poc_exact(model,q)
  analytical <- poc_bounds(model,q,use_zig=TRUE)
  constraints <- lapply(c("0","1"),function(g)
    poc_constraint(poc_query(setNames("1",g)),
                   max(0,risk[g]-.05),min(1,risk[g]+.05)))
  partial <- poc_partial_bounds(poc_partial_model(observed,
    constraints=constraints,lp_backend="zig"),q)
  # Cumulative paths only for the observed table; estimation uses annual y.
  cy <- t(apply(y,1,cummax))
  tab <- as.data.frame(table(history=as.character(group),
                            path=apply(cy,1,paste0,collapse="")))
  tab$prob <- tab$Freq/length(group)
  margins <- expand.grid(regime=c("0","1"),horizon=1:2,
                         stringsAsFactors=FALSE)
  margins$prob1 <- as.vector(cumulative$risks)
  trajectory <- poc_longitudinal_model(tab[,c("history","path","prob")],
    regimes=rbind("0"=c(0,0),"1"=c(1,1)), time_margins=margins,
    absorbing=TRUE,lp_backend="zig")
  curve <- poc_curve(trajectory,"1","0",bands=list(time_delta=.05))
  stopifnot(complete$status=="ok",partial$status=="ok",all(curve$status=="ok"),
    abs(complete$lower-analytical$lower)<1e-8,
    abs(complete$upper-analytical$upper)<1e-8,
    partial$lower<=complete$lower+1e-8,partial$upper>=complete$upper-1e-8)
  list(annual=annual,cumulative=cumulative,observed=observed,
       complete=complete,analytical=analytical,partial=partial,curve=curve)
}
pnsbib_adjusted_examples <- run_pnsbib_adjusted_examples()
