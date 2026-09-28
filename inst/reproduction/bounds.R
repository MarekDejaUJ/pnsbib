# Common-target fitted margins and factual path counts for 873 researchers.
source("common.R")
args<-commandArgs(TRUE);stopifnot(length(args)==2)
out<-init_output(args[1]);input<-args[2]
risks<-read.delim(file.path(input,"risks.tsv"))
radius<-read.delim(file.path(input,"bootstrap-summary.tsv"))$simultaneous_stability_radius
counts<-read.delim(file.path(input,"joint_counts.tsv"),colClasses=c(history="character",path="character"))
stopifnot(digest::digest(file=file.path(input,"joint_counts.tsv"),algo="sha256")==
 "d1bdd964be92d48eb2f50e4e66761ecea286ab9c07a45257b4e2375b61a8ada4")
regimes<-rbind("0"=rep(0,5),"1"=rep(1,5))
rows<-checks<-witnesses<-certificates<-list()
check<-function(label,ok) {
  checks[[length(checks)+1L]]<<-data.frame(check=label,passed=isTRUE(ok))
  if(!isTRUE(ok)) stop(label)
}
dir.create(file.path(out,"lp-evidence"))
.risk_certificate<-function(frame,fit) {
  if(get("lp_backend",frame)!="zig") return(invisible(NULL))
  A<-get("A",frame);b<-get("rhs",frame);c<-get("objective",frame)
  dirs<-match(get("directions",frame),c("<=","=",">="))-2L
  sense<-if(get("direction",frame)=="max") 1 else -1
  id<-length(certificates)+1L
  if(fit$status==0L) {
    z<-fit$solution;u<-fit$dual;uz<-sense*u;az<-drop(A%*%z);residual<-az-b
    violation<-ifelse(dirs==0,abs(residual),dirs*(-residual))
    primal<-max(c(0,-z/(1+abs(z)),violation/(1+abs(b)+drop(abs(A)%*%abs(z)))))
    dual<-max(c(0,dirs*uz/(1+abs(uz)),
      (sense*c-drop(crossprod(A,uz)))/(1+abs(c)+drop(crossprod(abs(A),abs(uz))))))
    gap<-max(abs(sum(c*z)-fit$objval),abs(sum(b*u)-fit$objval))/(1+abs(fit$objval))
    stopifnot(all(is.finite(c(primal,dual,gap))),max(primal,dual,gap)<=1e-8)
  } else primal<-dual<-gap<-NA_real_
  filename<-sprintf("lp-%05d.rds",id)
  saveRDS(list(A=A,b=b,c=c,dirs=dirs,sense=sense,fit=fit),file.path(out,"lp-evidence",filename),version=3,compress="gzip")
  certificates[[id]]<<-data.frame(id=id,status=fit$status,objective=fit$objval,
    rows=nrow(A),variables=ncol(A),primal_error=primal,dual_error=dual,gap=gap,
    file=filename,sha256=digest::digest(file=file.path(out,"lp-evidence",filename),algo="sha256"))
}
trace(".poc_lp",where=asNamespace("pnsbib"),exit=quote(.risk_certificate(environment(),returnValue())),print=FALSE)
record<-function(mode,h,kind,e,info,method,fit) {
  rows[[length(rows)+1L]]<<-data.frame(outcome=mode,year=2020+h,horizon=h,estimand=kind,
    epsilon=e,information=info,method=method,lower=fit$lower,upper=fit$upper,status=fit$status)
  check(paste("status",mode,h,kind,e,info,method),fit$status %in% c("ok","infeasible","undefined_condition"))
  if(fit$status=="ok") check(paste("probability_range",mode,h,kind,e,method),
    all(is.finite(ends(fit)))&&fit$lower>=-1e-8&&fit$upper<=1+1e-8&&fit$lower<=fit$upper+1e-8)
}
for(mode in c("annual","cumulative")) {
  rtab<-risks[risks$outcome==mode,];r<-sapply(2021:2025,function(y) rtab$adjusted[rtab$year==y][order(rtab$group[rtab$year==y])])
  rownames(r)<-c("0","1")
  z<-counts[counts$outcome==mode,];obs<-transform(z,prob=Freq/sum(Freq))
  tm<-expand.grid(regime=c("0","1"),horizon=1:5,stringsAsFactors=FALSE);tm$prob1<-as.vector(r)
  free<-poc_longitudinal_model(obs,regimes,absorbing=mode=="cumulative",lp_backend="zig")
  center<-poc_longitudinal_model(obs,regimes,time_margins=tm,absorbing=mode=="cumulative",lp_backend="zig")
  curves<-list(free=poc_curve(free,"1","0"))
  widths<-c(0,.025,.05,.1,.25,1,radius)
  for(k in seq_along(widths)) {
    e<-widths[k];key<-paste0("band",k)
    curves[[key]]<-poc_curve(center,"1","0",bands=list(time_delta=e))
    ref<-poc_curve(center,"1","0",bands=list(time_delta=e),lp_backend="reference")
    check(paste("native_reference_curve",mode,k),same(ends(curves[[key]]),ends(ref))&&identical(curves[[key]]$status,ref$status))
    cat(mode,"curve",k,"of",length(widths),"complete\n");flush.console()
  }
  for(h in 1:5) {
    n<-tapply(z$Freq,z$history,sum);success<-tapply(z$Freq*as.integer(substr(z$path,h,h)),z$history,sum)
    o<-cbind("0"=(n-success)/873,"1"=success/873)
    full<-poc_model(o,cbind("0"=1-r[,h],"1"=r[,h]),lp_backend="zig")
    for(kind in names(queries)) {
      q<-queries[[kind]]
      take<-function(key) curves[[key]][curves[[key]]$horizon==h & curves[[key]]$estimand==kind,,drop=FALSE]
      f0<-poc_partial_bounds(poc_partial_model(o,lp_backend="zig"),q)
      record(mode,h,kind,NA_real_,"consistency_only","partial_exact_lp",f0)
      record(mode,h,kind,NA_real_,"consistency_only","trajectory_exact_lp",take("free"))
      methods<-list(li_pearl=poc_bounds(full,q,use_zig=TRUE),shu=poc_shu2026(full,q,use_zig=TRUE),exact_lp=poc_exact(full,q))
      for(method in names(methods)) {
        record(mode,h,kind,0,"complete_adjusted_margins",method,methods[[method]])
        check(paste("complete_allocation",mode,h,kind,method),same(ends(methods[[method]]),binary_oracle(o,r[,h],0,kind)))
      }
      for(k in seq_along(widths)) {
        e<-widths[k];info<-if(k==length(widths)) "bootstrap_input_stability_box" else "population_risk_bands"
        pm<-poc_partial_model(o,risk_constraints(r[,h],e),lp_backend="zig")
        fit<-poc_partial_bounds(pm,q,return_witness=TRUE)
        ref<-poc_partial_bounds(pm,q,lp_backend="reference");joint<-take(paste0("band",k))
        record(mode,h,kind,e,info,"partial_exact_lp",fit)
        record(mode,h,kind,e,info,"trajectory_exact_lp",joint)
        check(paste("native_reference_static",mode,h,kind,k),same(ends(fit),ends(ref)))
        check(paste("allocation_static",mode,h,kind,k),same(ends(fit),binary_oracle(o,r[,h],e,kind)))
        check(paste("nested_joint",mode,h,kind,k),joint$lower>=fit$lower-1e-8&&joint$upper<=fit$upper+1e-8)
        for(side in c("lower","upper")) {
          w<-attr(fit,"witnesses")[[side]]
          fy<-ifelse(w$factual_x=="0",w$Y_0,w$Y_1)
          wo<-matrix(0,2,2,dimnames=dimnames(o))
          for(j in seq_len(nrow(w))) wo[w$factual_x[j],fy[j]]<-wo[w$factual_x[j],fy[j]]+w$prob[j]
          wr<-c(sum(w$prob[w$Y_0=="1"]),sum(w$prob[w$Y_1=="1"]))
          sel<-rep(TRUE,nrow(w));den<-1
          if(kind=="PN") {sel<-w$factual_x=="1"&fy=="1";den<-o[2,2]}
          if(kind=="PS") {sel<-w$factual_x=="0"&fy=="0";den<-o[1,1]}
          value<-sum(w$prob[sel&w$Y_0=="0"&w$Y_1=="1"])/den
          check(paste("endpoint",mode,h,kind,k,side),min(w$prob)>=-1e-8&&abs(sum(w$prob)-1)<1e-8&&
            max(abs(wo-o))<1e-8&&all(wr>=pmax(0,r[,h]-e)-1e-8)&&all(wr<=pmin(1,r[,h]+e)+1e-8)&&abs(value-fit[[side]])<1e-8)
          witnesses[[length(witnesses)+1L]]<-cbind(outcome=mode,year=2020+h,estimand=kind,epsilon=e,endpoint=side,w)
        }
      }
    }
  }
}
untrace(".poc_lp",where=asNamespace("pnsbib"))
save_tsv(do.call(rbind,rows),"bounds");save_tsv(do.call(rbind,checks),"checks")
save_tsv(do.call(rbind,witnesses),"endpoint-witnesses");save_tsv(do.call(rbind,certificates),"lp-certificates")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("All adjusted attribution bounds and witnesses verified\n")
