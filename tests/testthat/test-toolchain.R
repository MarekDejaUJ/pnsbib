toolchain <- new.env(parent=baseenv())
sys.source(system.file("toolchain","zig-toolchain.R",package="pnsbib"),envir=toolchain)

test_that("default compiler storage is session-local and persistent storage is opt-in", {
  expect_identical(toolchain$pnsbib_zig_cache(""),file.path(tempdir(),"pnsbib-zig"))
  expect_identical(toolchain$pnsbib_zig_cache("/explicit/cache"),"/explicit/cache/zig")
})

test_that("the compiler stays alive through build and temporary storage is cleaned", {
  old <- Sys.getenv("PNSBIB_ZIG_CACHE",unset=NA_character_)
  Sys.unsetenv("PNSBIB_ZIG_CACHE")
  on.exit(if(is.na(old)) Sys.unsetenv("PNSBIB_ZIG_CACHE") else
    Sys.setenv(PNSBIB_ZIG_CACHE=old),add=TRUE)
  cache <- toolchain$pnsbib_zig_cache("")
  resolver <- function() {
    dir.create(cache,showWarnings=FALSE)
    path <- file.path(cache,"zig");writeLines("fixture",path);path
  }
  runner <- function(command,args) {
    expect_true(file.exists(command))
    expect_identical(Sys.getenv("TMPDIR"),tempdir())
    expect_identical(args,c("build",shQuote(c("-j2","-Dname=with space"))))
    0L
  }
  expect_message(toolchain$pnsbib_zig_build(c("-j2","-Dname=with space"),runner,resolver,
                                            probe=function(path) "0.17.0"),"Compiling with Zig 0[.]17[.]0: ")
  expect_false(dir.exists(cache))
  expect_error(suppressMessages(toolchain$pnsbib_zig_build(character(),
    function(...) 1L,resolver,probe=function(path) "0.16.0")),"compilation failed")
  expect_false(dir.exists(cache))
})

test_that("automatic compiler archives have fixed platform-specific checksums", {
  hashes <- character()
  for(version in c("0.16.0","0.17.0")) for(os in c("Darwin","Linux","Windows"))
    for(arch in c("aarch64","x86_64")) {
      s <- toolchain$pnsbib_zig_spec(os,arch,version)
      expect_identical(s$version,version)
      expect_match(s$sha256,"^[a-f0-9]{64}$")
      expect_identical(s$url,paste0("https://ziglang.org/download/",version,"/zig-",s$key,"-",
                                    version,".",s$format))
      expect_gt(s$size,50000000)
      expect_identical(s$format,if(os=="Windows") "zip" else "tar.xz")
      hashes <- c(hashes,s$sha256)
    }
  expect_identical(anyDuplicated(hashes),0L)
  expect_identical(toolchain$pnsbib_zig_spec("Darwin","arm64","0.16.0")$sha256,
                   "b23d70deaa879b5c2d486ed3316f7eaa53e84acf6fc9cc747de152450d401489")
  expect_identical(toolchain$pnsbib_zig_spec("Darwin","arm64","0.17.0")$sha256,
                   "b607e9b9234790a008116ae5bdb71c6243b84b9fb42a53a9e70fde41c06c536a")
  expect_error(toolchain$pnsbib_zig_spec("Linux","x86_64","0.15.2"),"No pinned Zig archive")
  expect_identical(toolchain$pnsbib_zig_spec("Darwin","arm64"),
                   toolchain$pnsbib_zig_spec("Darwin","aarch64"))
  expect_identical(toolchain$pnsbib_zig_spec("Linux","amd64"),
                   toolchain$pnsbib_zig_spec("Linux","x86_64"))
  expect_error(toolchain$pnsbib_zig_spec("Plan9","x86_64"),"supports macOS")
  expect_error(toolchain$pnsbib_zig_spec("Linux","armv7"),"supports aarch64")
})

test_that("Zig 0.16 and 0.17 are accepted and other series are rejected", {
  for(v in c("0.16.0","0.16.1","0.17.0","0.17.2"))
    expect_true(toolchain$pnsbib_zig_supported(v))
  for(v in c("","0.15.2","0.18.0","0.16.0-dev.1+abc","1.16.0","0.160.0","zig 0.16.0"))
    expect_false(toolchain$pnsbib_zig_supported(v))
})

test_that("automatic setup downloads Zig 0.16.0 unless 0.17.0 is requested", {
  expect_identical(toolchain$pnsbib_zig_download_version(""),"0.16.0")
  expect_identical(toolchain$pnsbib_zig_download_version("0.16.0"),"0.16.0")
  expect_identical(toolchain$pnsbib_zig_download_version("0.17.0"),"0.17.0")
  expect_error(toolchain$pnsbib_zig_download_version("0.18.0"),"0.16.0 or 0.17.0; found 0.18.0")
})

test_that("an explicit ZIG executable must belong to a supported series", {
  path <- tempfile("zig-fixture-");writeLines("fixture",path)
  old <- Sys.getenv("ZIG",unset=NA_character_)
  probe <- toolchain$pnsbib_zig_version
  on.exit({unlink(path);assign("pnsbib_zig_version",probe,envir=toolchain)
    if(is.na(old)) Sys.unsetenv("ZIG") else Sys.setenv(ZIG=old)},add=TRUE)
  Sys.setenv(ZIG=path)
  for(v in c("0.16.0","0.17.0")) {
    toolchain$pnsbib_zig_version <- function(path) v
    expect_identical(toolchain$pnsbib_zig_resolve(),normalizePath(path,winslash="/"))
  }
  for(v in c("0.15.2","0.18.0")) {
    toolchain$pnsbib_zig_version <- function(path) v
    expect_error(toolchain$pnsbib_zig_resolve(),paste0("found ",v))
  }
  toolchain$pnsbib_zig_version <- function(path) ""
  expect_error(toolchain$pnsbib_zig_resolve(),"found no working Zig")
})

test_that("archive paths cannot escape the verified root", {
  expect_true(toolchain$pnsbib_zig_members(c("zig-test/","zig-test/zig","zig-test/lib/std.zig"),"zig-test"))
  for(path in c("../zig","zig-test/../escape","/zig-test/zig","C:/zig-test/zig",
               "another-root/zig","zig-test\\zig","zig-test/lib/../../escape"))
    expect_error(toolchain$pnsbib_zig_members(path,"zig-test"),"Unsafe")
  expect_error(toolchain$pnsbib_zig_members(character(),"zig-test"),"Unsafe")
  expect_error(toolchain$pnsbib_zig_members(NA_character_,"zig-test"),"Unsafe")
})

test_that("version probes leave command quoting to system2", {
  path <- tempfile("compiler path ");writeLines("fixture",path)
  on.exit(unlink(path))
  toolchain$system2 <- function(command,args,...) {
    expect_identical(command,path)
    expect_identical(args,"version")
    "0.16.0"
  }
  on.exit(rm("system2",envir=toolchain),add=TRUE)
  expect_identical(toolchain$pnsbib_zig_version(path),"0.16.0")
  expect_identical(toolchain$pnsbib_zig_version(paste0(path,"-missing")),"")
})

toolchain_fixture <- function() {
  root <- tempfile("toolchain-test-");dir.create(root)
  top <- "zig-test-0.16.0";dir.create(file.path(root,top))
  writeLines("a compiler fixture, never executed",file.path(root,top,"zig"))
  archive <- file.path(root,"fixture.tar")
  previous <- setwd(root);on.exit(setwd(previous))
  utils::tar(archive,files=top,compression="none",tar="internal")
  list(root=root,archive=archive,cache=file.path(root,"cache"),
       spec=list(version="0.16.0",key="test",top=top,format="tar.xz",executable="zig",
                 url="https://ziglang.org/download/fixture",size=file.info(archive)$size,
                 sha256=digest::digest(file=archive,algo="sha256")))
}

test_that("offline extraction and cache reuse need neither network nor a compiler on PATH", {
  f <- toolchain_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  never_download <- function(...) stop("unexpected download")
  probe <- function(path) "0.16.0"
  expect_message(path <- toolchain$pnsbib_zig_install(f$spec,f$cache,f$archive,FALSE,
                                                     never_download,probe),"Verified Zig")
  expect_true(file.exists(path))
  expect_true(file.exists(file.path(dirname(path),".pnsbib-toolchain.dcf")))
  expect_message(reused <- toolchain$pnsbib_zig_install(f$spec,f$cache,download=FALSE,
                                                       downloader=never_download,probe=probe),"Reusing")
  expect_identical(path,reused)
  writeLines("changed binary",path)
  expect_error(toolchain$pnsbib_zig_install(f$spec,f$cache,probe=probe),"integrity check failed")
})

test_that("downloaded bytes are checked before extraction or execution", {
  f <- toolchain_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  calls <- 0L
  downloader <- function(url,destfile,...) {
    expect_identical(url,f$spec$url);calls <<- calls+1L
    stopifnot(file.copy(f$archive,destfile));0L
  }
  expect_message(path <- toolchain$pnsbib_zig_install(f$spec,f$cache,downloader=downloader,
                                                     probe=function(path) "0.16.0"),"Downloading")
  expect_identical(calls,1L)
  expect_true(file.exists(path))
  bad <- f$spec;bad$sha256 <- paste(rep("0",64),collapse="")
  expect_error(toolchain$pnsbib_zig_install(bad,file.path(f$root,"bad"),f$archive,
                     probe=function(path) stop("must not execute")),"SHA-256 mismatch")
  expect_false(dir.exists(file.path(f$root,"bad",bad$top)))
  expect_false(dir.exists(file.path(f$root,"bad",paste0(bad$top,".lock"))))
  bad$size <- bad$size+1
  expect_error(toolchain$pnsbib_zig_install(bad,file.path(f$root,"size"),f$archive),"size or SHA-256")
})

test_that("disabled downloads, incomplete caches and concurrent installation fail clearly", {
  f <- toolchain_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  expect_error(toolchain$pnsbib_zig_install(f$spec,f$cache,download=FALSE),"download is disabled")
  expect_false(dir.exists(f$cache))
  expect_error(toolchain$pnsbib_zig_install(f$spec,"relative-cache",f$archive),"absolute path")
  expect_error(toolchain$pnsbib_zig_install(f$spec,file.path(getwd(),"toolchain-cache"),f$archive),"outside")
  dir.create(f$cache)
  dir.create(paste0(file.path(f$cache,f$spec$top),".lock"))
  expect_error(toolchain$pnsbib_zig_install(f$spec,f$cache,f$archive),"locked")
  dir.create(file.path(f$cache,f$spec$top))
  expect_error(toolchain$pnsbib_zig_install(f$spec,f$cache,f$archive),"Incomplete Zig cache")
  expect_error(toolchain$pnsbib_zig_install(f$spec,file.path(f$root,"network"),
               downloader=function(...) stop("network unavailable")),"Zig download failed")
  expect_error(toolchain$pnsbib_zig_install(f$spec,file.path(f$root,"status"),
               downloader=function(...) 1L),"did not complete")
  expect_error(toolchain$pnsbib_zig_install(f$spec,file.path(f$root,"version"),f$archive,
               probe=function(path) "0.15.2"),"required Zig executable")
})
