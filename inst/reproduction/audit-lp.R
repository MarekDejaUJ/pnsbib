# Independent verification of saved native optimization certificates, no solver calls.
args<-commandArgs(TRUE);stopifnot(length(args)==2)
input<-args[1];out<-args[2];stopifnot(!file.exists(out),dir.exists(dirname(out)))
dir.create(out,recursive=TRUE)
ledger<-read.delim(file.path(input,"lp-certificates.tsv"));rows<-list()
for(i in seq_len(nrow(ledger))) {
  file<-file.path(input,"lp-evidence",ledger$file[i])
  stopifnot(digest::digest(file=file,algo="sha256")==ledger$sha256[i])
  z<-readRDS(file);f<-z$fit
  stopifnot(f$status==0L)
  ax<-drop(z$A%*%f$solution);atv<-drop(crossprod(z$A,z$sense*f$dual))
  rel_primal<-max(c(0,-f$solution/(1+abs(f$solution)),
    ifelse(z$dirs==0,abs(ax-z$b),z$dirs*(z$b-ax))/(1+abs(z$b)+drop(abs(z$A)%*%abs(f$solution)))))
  rel_dual<-max(c(0,z$dirs*z$sense*f$dual/(1+abs(f$dual)),
    (z$sense*z$c-atv)/(1+abs(z$c)+drop(crossprod(abs(z$A),abs(f$dual))))))
  value<-sum(z$c*f$solution);dualvalue<-sum(z$b*f$dual)
  gap<-max(abs(value-f$objval),abs(dualvalue-f$objval))/(1+abs(f$objval))
  stopifnot(all(is.finite(c(rel_primal,rel_dual,gap))),max(rel_primal,rel_dual,gap)<=1e-8)
  rows[[i]]<-data.frame(id=i,primal_error=rel_primal,dual_error=rel_dual,gap=gap,passed=TRUE)
}
write.table(do.call(rbind,rows),file.path(out,"independent-lp-audit.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat(nrow(ledger),"native programs verified independently\n")
