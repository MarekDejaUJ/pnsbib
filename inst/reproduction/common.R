# Fixed-design numerical validation and independent binary allocation oracle.
library(pnsbib)
stopifnot(as.character(packageVersion("pnsbib"))=="0.2.0")
base_dir <- "."
save_tsv <- function(x,name) write.table(x,file.path(out,paste0(name,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
queries <- list(PNS=poc_query(c("0"="0","1"="1")),
  PN=poc_query(c("0"="0"),"1","1",TRUE),
  PS=poc_query(c("1"="1"),"0","0",TRUE))
ends <- function(x) c(x$lower,x$upper)
same <- function(a,b) isTRUE(all.equal(a,b,tolerance=1e-8,check.attributes=FALSE))
observed_cells <- function(g,y) {
  o <- unclass(table(factor(g,0:1),factor(y,0:1)))/length(g)
  dimnames(o) <- list(c("0","1"),c("0","1")); o
}
scaled_design <- function(h,target=h) {
  center <- colMeans(h); scale <- apply(h,2,sd); scale[!is.finite(scale)|scale==0] <- 1
  list(x=cbind(intercept=1,sweep(sweep(h,2,center),2,scale,"/")),
       target=cbind(intercept=1,sweep(sweep(target,2,center),2,scale,"/")),center=center,scale=scale)
}
risk_constraints <- function(r,e) lapply(1:2,function(j) poc_constraint(
  poc_query(setNames("1",as.character(j-1))), max(0,r[j]-e),min(1,r[j]+e)))
# Independent two-binary-variable allocation oracle; no package internals.
binary_oracle <- function(o,r,e,kind) {
  p <- rowSums(o)
  a <- c(max(0,r[2]-e-o[2,2]),min(p[1],r[2]+e-o[2,2]))
  b <- c(max(0,r[1]-e-o[1,2]),min(p[2],r[1]+e-o[1,2]))
  if (a[1]>a[2]+1e-9 || b[1]>b[2]+1e-9) return(c(NA_real_,NA_real_))
  active <- c(max(0,o[2,2]-b[2]),min(o[2,2],p[2]-b[1]))
  reference <- c(max(0,a[1]-o[1,2]),min(o[1,1],a[2]))
  if (kind=="PN" && o[2,2]==0 || kind=="PS" && o[1,1]==0) return(c(NA_real_,NA_real_))
  unname(switch(kind,PNS=active+reference,PN=active/o[2,2],PS=reference/o[1,1]))
}
init_output <- function(path) {
  stopifnot(!file.exists(path),dir.exists(dirname(path)))
  dir.create(path,recursive=TRUE);normalizePath(path)
}
