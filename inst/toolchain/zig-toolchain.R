# Source-installation helper. Compiler archives are pinned to the publisher's
# SHA-256 values at https://ziglang.org/download/index.json (Zig 0.16.0).
pnsbib_zig_spec <- function(os=Sys.info()[["sysname"]], arch=R.version$arch) {
  os <- switch(tolower(os),darwin="macos",linux="linux",windows="windows",
               stop("Automatic Zig setup supports macOS, Linux and Windows."))
  arch <- switch(tolower(arch),arm64="aarch64",aarch64="aarch64",
                 x86_64="x86_64",amd64="x86_64",
                 stop("Automatic Zig setup supports aarch64 and x86_64 R."))
  key <- paste(arch,os,sep="-")
  hashes <- c(
    "aarch64-macos"="b23d70deaa879b5c2d486ed3316f7eaa53e84acf6fc9cc747de152450d401489",
    "x86_64-macos"="0387557ed1877bc6a2e1802c8391953baddba76081876301c522f52977b52ba7",
    "aarch64-linux"="ea4b09bfb22ec6f6c6ceac57ab63efb6b46e17ab08d21f69f3a48b38e1534f17",
    "x86_64-linux"="70e49664a74374b48b51e6f3fdfbf437f6395d42509050588bd49abe52ba3d00",
    "aarch64-windows"="aee38316ee4111717900f45dd3130145c39289e105541d737eb8c5ed653c78ef",
    "x86_64-windows"="68659eb5f1e4eb1437a722f1dd889c5a322c9954607f5edcf337bc3684a75a7e")
  sizes <- c("aarch64-macos"=52238004,"x86_64-macos"=57396836,
             "aarch64-linux"=51211944,"x86_64-linux"=55478392,
             "aarch64-windows"=93109828,"x86_64-windows"=97217739)
  top <- paste0("zig-",key,"-0.16.0")
  format <- if(os=="windows") "zip" else "tar.xz"
  list(version="0.16.0",key=key,top=top,format=format,
       executable=if(os=="windows") "zig.exe" else "zig",
       url=paste0("https://ziglang.org/download/0.16.0/",top,".",format),
       sha256=unname(hashes[key]),size=unname(sizes[key]))
}

pnsbib_zig_version <- function(path) {
  if(!file.exists(path) || dir.exists(path)) return("")
  answer <- tryCatch(suppressWarnings(system2(path,"version",stdout=TRUE,
                                             stderr=TRUE,timeout=15)),error=function(e) "")
  if(!is.null(attr(answer,"status")) || length(answer)!=1L) return("")
  trimws(answer)
}

pnsbib_zig_members <- function(members,top) {
  if(!length(members) || anyNA(members) || any(grepl("\\\\",members)) ||
     any(grepl("(^|/)\\.\\.(/|$)|^/|^[A-Za-z]:",members)) ||
     any(!(members==top | startsWith(members,paste0(top,"/")))))
    stop("Unsafe or unexpected paths in the Zig archive.")
  invisible(TRUE)
}

pnsbib_zig_install <- function(spec,cache,archive="",download=TRUE,
                               downloader=utils::download.file,probe=pnsbib_zig_version) {
  if(!requireNamespace("digest",quietly=TRUE))
    stop("The digest R dependency is required to verify the Zig archive.")
  cache <- path.expand(cache)
  if(!grepl("^(/|[A-Za-z]:[/\\\\])",cache)) stop("The Zig cache must be an absolute path.")
  # Resolve existing ancestors before creating anything; never build in the source tree.
  ancestor <- cache; suffix <- character()
  while(!file.exists(ancestor)) {
    suffix <- c(basename(ancestor),suffix); ancestor <- dirname(ancestor)
  }
  cache <- do.call(file.path,as.list(c(normalizePath(ancestor,winslash="/",mustWork=TRUE),suffix)))
  source <- normalizePath(getwd(),winslash="/",mustWork=TRUE)
  if(identical(cache,source) || startsWith(cache,paste0(source,"/")))
    stop("The Zig cache must be outside the package source directory.")
  target <- file.path(cache,spec$top)
  binary <- file.path(target,spec$executable)
  receipt <- file.path(target,".pnsbib-toolchain.dcf")
  if(dir.exists(target)) {
    if(!file.exists(receipt) || !file.exists(binary))
      stop("Incomplete Zig cache entry: ",target,". Use a fresh PNSBIB_ZIG_CACHE directory.")
    r <- read.dcf(receipt)
    if(!all(c("Archive-SHA256","Binary-SHA256") %in% colnames(r)) ||
       r[1,"Archive-SHA256"]!=spec$sha256 ||
       r[1,"Binary-SHA256"]!=digest::digest(file=binary,algo="sha256"))
      stop("Zig cache integrity check failed: ",target)
    if(probe(binary)!=spec$version) stop("Cached Zig version does not match ",spec$version)
    message("Reusing cached Zig ",spec$version)
    return(normalizePath(binary,winslash="/",mustWork=TRUE))
  }
  if(!nzchar(archive) && !download)
    stop("No Zig ",spec$version," found; automatic download is disabled. Supply ZIG or PNSBIB_ZIG_ARCHIVE.")
  if(!dir.exists(cache) && !dir.create(cache,recursive=TRUE)) stop("Cannot create Zig cache: ",cache)
  lock <- paste0(target,".lock")
  if(!dir.create(lock,showWarnings=FALSE))
    stop("Zig cache is locked by another installation: ",lock,". Retry when it finishes.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  stage <- tempfile("pnsbib-zig-",tmpdir=cache)
  if(!dir.create(stage)) stop("Cannot create toolchain staging directory.")
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  if(!nzchar(archive)) {
    archive <- file.path(stage,paste0(spec$top,".",spec$format))
    message("Downloading Zig ",spec$version," to the user cache (",round(spec$size/1024^2)," MiB).")
    old <- options(timeout=max(300,getOption("timeout",60)))
    on.exit(options(old),add=TRUE)
    status <- tryCatch(downloader(spec$url,archive,mode="wb",quiet=TRUE,method="libcurl"),
      error=function(e) stop("Zig download failed. Check network access or supply PNSBIB_ZIG_ARCHIVE. ",conditionMessage(e)))
    if(!identical(as.integer(status),0L)) stop("Zig download did not complete.")
  }
  archive <- normalizePath(archive,winslash="/",mustWork=TRUE)
  if(file.info(archive)$size!=spec$size || digest::digest(file=archive,algo="sha256")!=spec$sha256)
    stop("Zig archive size or SHA-256 mismatch; nothing was extracted or executed.")
  unpack <- file.path(stage,"unpacked");dir.create(unpack)
  if(spec$format=="zip") {
    members <- utils::unzip(archive,list=TRUE)$Name
    pnsbib_zig_members(members,spec$top)
    utils::unzip(archive,exdir=unpack)
  } else {
    members <- utils::untar(archive,list=TRUE)
    pnsbib_zig_members(members,spec$top)
    status <- utils::untar(archive,exdir=unpack)
    if(!identical(as.integer(status),0L)) stop("Could not extract the verified Zig archive.")
  }
  extracted <- file.path(unpack,spec$top)
  executable <- file.path(extracted,spec$executable)
  if(!file.exists(executable) || nzchar(Sys.readlink(executable)) || probe(executable)!=spec$version)
    stop("The verified archive did not provide the required Zig executable.")
  receipt_data <- matrix(c(spec$version,spec$url,spec$sha256,digest::digest(file=executable,algo="sha256")),
                         nrow=1,dimnames=list(NULL,c("Version","URL","Archive-SHA256","Binary-SHA256")))
  write.dcf(receipt_data,file.path(extracted,".pnsbib-toolchain.dcf"))
  if(!file.rename(extracted,target)) stop("Could not finalize Zig cache: ",target)
  message("Verified Zig ",spec$version," is ready in the user cache.")
  normalizePath(binary,winslash="/",mustWork=TRUE)
}

pnsbib_zig_resolve <- function() {
  version <- "0.16.0"
  explicit <- Sys.getenv("ZIG","")
  if(nzchar(explicit)) {
    if(!file.exists(explicit)) explicit <- unname(Sys.which(explicit))
    if(!nzchar(explicit) || pnsbib_zig_version(explicit)!=version)
      stop("ZIG must name a working Zig ",version," executable; unset ZIG for automatic setup.")
    return(normalizePath(explicit,winslash="/",mustWork=TRUE))
  }
  found <- unname(Sys.which("zig"))
  if(nzchar(found) && pnsbib_zig_version(found)==version)
    return(normalizePath(found,winslash="/",mustWork=TRUE))
  spec <- pnsbib_zig_spec()
  cache <- file.path(Sys.getenv("PNSBIB_ZIG_CACHE",tools::R_user_dir("pnsbib","cache")),"zig")
  archive <- Sys.getenv("PNSBIB_ZIG_ARCHIVE","")
  download <- !tolower(Sys.getenv("PNSBIB_ZIG_DOWNLOAD","true")) %in% c("0","false","no")
  pnsbib_zig_install(spec,cache,archive,download)
}

if(sys.nframe()==0L) {
  tryCatch(cat(pnsbib_zig_resolve(),"\n",sep=""),error=function(e) {
    message(conditionMessage(e));quit(status=1L)
  })
}
