# Plot verified attribution intervals; no estimation or optimization is performed.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L, file.exists(args[1]), !file.exists(args[2]),
          dir.exists(dirname(args[2])))
b <- read.delim(args[1])
stopifnot(nrow(b)==570L, all(b$status=="ok"), all(is.finite(b$lower)),
          all(is.finite(b$upper)), all(b$lower>=-1e-8), all(b$upper<=1+1e-8),
          all(b$lower<=b$upper+1e-8))
selected <- b$information=="population_risk_bands" &
  !is.na(b$epsilon) & b$epsilon %in% c(0,.05)
z <- b[selected & b$method=="partial_exact_lp",]
key <- function(x) paste(x$outcome,x$year,x$estimand,x$epsilon,sep=":")
stopifnot(nrow(z)==60L, !anyDuplicated(key(z)))
j <- b[selected & b$method=="trajectory_exact_lp",]
stopifnot(nrow(j)==60L, setequal(key(z),key(j)), !anyDuplicated(key(j)))
j <- j[match(key(z),key(j)),]
delta <- max(abs(as.matrix(z[c("lower","upper")])-as.matrix(j[c("lower","upper")])))
stopifnot(delta<1e-8)
exact <- z[z$epsilon==0,]
for(method in c("li_pearl","shu","exact_lp")) {
  m <- b[b$information=="complete_adjusted_margins" & b$method==method,]
  stopifnot(nrow(m)==30L, setequal(key(exact),key(m)))
  m <- m[match(key(exact),key(m)),]
  stopifnot(max(abs(as.matrix(exact[c("lower","upper")])-
                       as.matrix(m[c("lower","upper")])))<1e-8)
}
band <- z[z$epsilon==.05,]
band$epsilon <- 0
band <- band[match(key(exact),key(band)),]
stopifnot(all(band$lower<=exact$lower+1e-8), all(band$upper>=exact$upper-1e-8))
out <- args[2]
dir.create(out)
write.table(z,file.path(out,"plot-data.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(plotted_intervals=nrow(z),static_trajectory_max_difference=delta,
                       complete_exact_agreement=TRUE,nested_bands=TRUE),
            file.path(out,"plot-checks.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
draw <- function(gray=FALSE) {
  layout(rbind(1:2,3:4,5:6,c(7,7)),heights=c(1,1,1,.34))
  par(mar=c(3.2,3.5,1.8,.6),oma=c(0,0,.2,0),mgp=c(2.15,.6,0),
      family="sans",cex=.90,las=1,tcl=-.2)
  cols <- if(gray) c("#1A1A1A","#777777") else c("#0072B2","#D55E00")
  k <- 0L
  for(kind in c("PNS","PN","PS")) for(outcome in c("annual","cumulative")) {
    k <- k+1L
    plot(NA,xlim=c(2020.7,2025.3),ylim=c(-.025,1.025),xaxt="n",yaxt="n",
         xlab=if(kind=="PS") "Year" else "",ylab=paste(kind,"probability"),
         bty="l",yaxs="i")
    axis(1,at=2021:2025,labels=paste0("'",21:25))
    axis(2,at=seq(0,1,.25))
    abline(h=seq(0,1,.25),col="#E5E5E5",lwd=.5)
    for(a in 1:2) {
      v <- z[z$outcome==outcome & z$estimand==kind & z$epsilon==c(0,.05)[a],]
      v <- v[order(v$year),]
      stopifnot(nrow(v)==5L, identical(as.integer(v$year),2021:2025))
      xx <- v$year+c(-.08,.08)[a]
      segments(xx,v$lower,xx,v$upper,col=cols[a],lty=a,lwd=1)
      for(side in c("lower","upper")) {
        lines(xx,v[[side]],col=cols[a],lty=a,lwd=.7)
        points(xx,v[[side]],col=cols[a],pch=c(16,1)[a],cex=.55)
      }
    }
    title(main=paste(LETTERS[k],if(outcome=="annual") "Annual" else "Cumulative"),
          cex.main=1,font.main=2)
  }
  par(mar=c(0,0,0,0));plot.new()
  legend(.5,.9,xjust=.5,yjust=1,
         legend=c("Exact adjusted risks","Adjusted risk bands +/- 0.05"),
         col=cols,pch=c(16,1),lty=1:2,bty="n",cex=.92,y.intersp=1.1)
}
if(capabilities("aqua")) {
  quartz(type="pdf",file=file.path(out,"curves.pdf"),width=117/25.4,height=6.1,pointsize=10)
} else if(capabilities("cairo")) {
  cairo_pdf(file.path(out,"curves.pdf"),width=117/25.4,height=6.1,pointsize=10)
} else {
  pdf(file.path(out,"curves.pdf"),width=117/25.4,height=6.1,pointsize=10,useDingbats=FALSE)
}
draw();invisible(dev.off())
stopifnot(file.exists(file.path(out,"curves.pdf")), file.info(file.path(out,"curves.pdf"))$size>0)
png_type <- if(capabilities("aqua")) "quartz" else if(capabilities("cairo")) "cairo" else "Xlib"
png(file.path(out,"curves.png"),width=117/25.4,height=6.1,units="in",res=600,
    type=png_type,pointsize=10)
draw();invisible(dev.off())
png(file.path(out,"curves-grayscale.png"),width=117/25.4,height=6.1,units="in",res=150,
    type=png_type,pointsize=10)
draw(TRUE);invisible(dev.off())
cat("Verified and plotted 60 adjusted intervals in six panels.\n")
