# 600 repetitions, seed 2026092803; three mechanisms, 873 units per repetition.
source("common.R")
args <- commandArgs(TRUE);stopifnot(length(args)==1)
out <- init_output(args[1]);set.seed(2026092803)
risks <- bounds <- checks <- list();n<-873L
for (scenario in c("correct","misspecified","poor_overlap")) {
  for (rep in 1:200) {
    h <- runif(n,-2,2); z <- rbinom(n,1,.5)
    eta <- -.3+h+.3*z
    p0 <- plogis(-.6+1.2*h+.4*z);p1<-plogis(.2+1.2*h+.4*z)
    if(scenario=="misspecified") {
      eta <- -.3+.3*h+1.8*(h*h-4/3)+.3*z
      p0<-plogis(-.6+.3*h+2*(h*h-4/3)+.4*z)
      p1<-plogis(.2+.3*h+2*(h*h-4/3)+.4*z)
    }
    if(scenario=="poor_overlap") eta<-4*h+2*z-1
    e<-plogis(eta);g<-rbinom(n,1,e)
    y<-rbinom(n,1,ifelse(g==1,p1,p0))
    design<-scaled_design(cbind(h=h,z=z))
    fit<-poc_adjusted_risks(design$x,y,g)
    truth_risk<-c(mean(p0),mean(p1))
    crude<-as.numeric(tapply(y,g,mean));adjusted<-as.numeric(fit$risks)
    true_o<-matrix(c(mean((1-e)*(1-p0)),mean((1-e)*p0),mean(e*(1-p1)),mean(e*p1)),2,byrow=TRUE,
      dimnames=list(c("0","1"),c("0","1")))
    truth<-c(PNS=mean((1-p0)*p1),PN=mean(e*(1-p0)*p1)/mean(e*p1),
      PS=mean((1-e)*(1-p0)*p1)/mean((1-e)*(1-p0)))
    o<-observed_cells(g,y)
    for(method in c("crude","adjusted")) {
      r<-if(method=="crude") crude else adjusted
      for(j in 1:2) risks[[length(risks)+1L]]<-data.frame(scenario=scenario,replication=rep,method=method,
        group=j-1,estimate=r[j],truth=truth_risk[j],error=r[j]-truth_risk[j],
        poor_overlap_fraction=mean(e<.05|e>.95))
      for(width in c(0,.05)) for(kind in names(queries)) {
        b<-binary_oracle(o,r,width,kind);feasible<-all(is.finite(b))
        bounds[[length(bounds)+1L]]<-data.frame(scenario=scenario,replication=rep,method=method,
          width=width,estimand=kind,lower=b[1],upper=b[2],truth=truth[kind],feasible=feasible,
          contains=feasible && b[1]<=truth[kind]+1e-8 && b[2]>=truth[kind]-1e-8)
        # Every tenth replication also checks the package LP independently.
        if(rep %% 10==0) {
          model<-poc_partial_model(o,risk_constraints(r,width),lp_backend="zig")
          result<-poc_partial_bounds(model,queries[[kind]])
          ok<-if(feasible) result$status=="ok"&&same(ends(result),b) else result$status=="infeasible"
          checks[[length(checks)+1L]]<-data.frame(scenario=scenario,replication=rep,method=method,width=width,
            estimand=kind,check="native_vs_allocation",passed=ok)
          stopifnot(ok)
        }
      }
    }
    for(kind in names(queries)) {
      b<-binary_oracle(true_o,truth_risk,0,kind)
      ok<-b[1]<=truth[kind]+1e-10 && b[2]>=truth[kind]-1e-10
      checks[[length(checks)+1L]]<-data.frame(scenario=scenario,replication=rep,method="known_population",width=0,
        estimand=kind,check="population_truth_containment",passed=ok)
      stopifnot(ok)
    }
  }
  cat(scenario,"200 repetitions complete\n");flush.console()
}
save_tsv(do.call(rbind,risks),"simulation-risks")
save_tsv(do.call(rbind,bounds),"simulation-bounds")
save_tsv(do.call(rbind,checks),"simulation-checks")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("Simulation complete; all numerical truth/oracle checks passed\n")
