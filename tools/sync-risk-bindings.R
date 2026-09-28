# Restore internal native registrations after rzig::document().
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
pkg <- normalizePath(args[1],mustWork=TRUE)
manifest <- file.path(pkg,"src","rzig","framework","generated","manifest.zig")
namespace <- file.path(pkg,"NAMESPACE")
text <- paste(readLines(manifest,warn=FALSE),collapse="\n")
ns <- paste(readLines(namespace,warn=FALSE),collapse="\n")
entries <- list(
  poc_risk_fit_zig=list(doc="Internal penalized binary-risk fit.",
    parameters=c("n","p","m","x","y","weights","target","penalty","lambda","max_iterations","tolerance")),
  poc_risk_standardize_zig=list(doc="Internal target standardization.",
    parameters=c("n","t","pred","weights","cumulative")))
insert_after <- function(text,anchor,addition) {
  at <- gregexpr(anchor,text,fixed=TRUE)[[1]]
  stopifnot(length(at)==1L,at>0)
  end <- at+nchar(anchor)-1L
  paste0(substr(text,1,end),addition,substring(text,end+1L))
}
for(name in rev(names(entries))) {
  z <- entries[[name]]
  if(!grepl(paste0('.name = "',name,'"'),text,fixed=TRUE)) {
    entry <- paste0('\n            .{\n                .name = "',name,'",\n',
      '                .func = bound_root.',name,',\n',
      '                .doc = "',z$doc,'",\n',
      '                .parameters = .{ ',paste(paste0('"',z$parameters,'"'),collapse=', '),' },\n            },')
    text <- insert_after(text,'pub const exports = .{',entry)
  }
  binding <- paste0('  ',name,'_ = "',name,'",')
  if(!grepl(binding,ns,fixed=TRUE)) ns <- insert_after(ns,'useDynLib("pnsbib",',paste0('\n',binding))
}
writeLines(text,manifest,useBytes=TRUE)
writeLines(ns,namespace,useBytes=TRUE)
