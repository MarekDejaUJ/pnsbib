confidence_fixture <- function() {
  o <- matrix(c(.3,.2,.1,.4),2,byrow=TRUE,dimnames=list(c("0","1"),c("0","1")))
  a <- matrix(c(.7,.3,.4,.6),2,byrow=TRUE,dimnames=dimnames(o))
  model <- poc_model(o,a)
  p <- c(as.vector(t(o)),as.vector(t(a)))
  samples <- setNames(lapply(p,function(z) rep(c(0,1),c(100-round(100*z),round(100*z)))),
                      poc_margin_map(model)$cell_id)
  list(model=model,samples=samples)
}

confidence_make <- function(template, samples, method="hoeffding", ...)
  poc_confidence_region(template,samples,target="Synthetic prespecified population",
    sampling_unit="individual",method=method,assume_sampling=TRUE,...)

test_that("confidence cells match both explicit model systems", {
  f <- confidence_fixture()
  m <- poc_margin_map(f$model)
  expect_identical(m$type,rep(c("observed","intervention"),each=4))
  expect_identical(m$treatment,rep(rep(c("0","1"),each=2),2))
  expect_identical(m$outcome,rep(c("0","1"),4))
  o <- expand.grid(history=c("a","b"),path=c("ll","lh","hh"),stringsAsFactors=FALSE)
  o$prob <- 1/6
  regimes <- rbind(a=c(0,0),b=c(1,1))
  paths <- rbind(ll=c("low","low"),lh=c("low","high"),hh=c("high","high"))
  times <- data.frame(regime=c("b","a"),horizon=c(2,1),outcome=c("high","low"),prob=2/3)
  model <- poc_longitudinal_model(o,regimes,time_margins=times,outcome_paths=paths,
    absorbing=TRUE,outcome_order=c("low","high"))
  map <- poc_margin_map(model)
  expect_equal(nrow(map),8)
  expect_identical(map$history[1:6],rep(c("a","b"),each=3))
  expect_identical(map$path[1:6],rep(c("ll","lh","hh"),2))
  expect_identical(map$outcome[7:8],c("high","low"))
  r <- confidence_make(model,list())
  z <- poc_confidence_bounds(r,poc_longitudinal_query("pns","b","a",2,
    regime_outcome="high",reference_outcome="low"))
  expect_equal(c(z$lower,z$upper),c(0,1))
  expect_identical(z$outcome_order,"low;high")
  samples <- setNames(lapply(c(rep(1/6,6),2/3,2/3),function(p)
    rep(c(0,1),c(600-round(600*p),round(600*p)))),map$cell_id)
  r <- confidence_make(model,samples)
  delta <- sqrt(log(16/.05)/1200)
  for (kind in c("pns","pn","ps","persistent","event_time","lagged")) {
    q <- poc_longitudinal_query(kind,"b","a",2,regime_outcome="high",reference_outcome="low")
    actual <- poc_confidence_bounds(r,q)
    bands <- poc_longitudinal_sensitivity(model,q,observed_delta=delta,time_delta=delta)
    expect_equal(c(actual$lower,actual$upper),c(bands$lower,bands$upper),tolerance=1e-8)
    expect_identical(actual$status,bands$status)
  }
})

test_that("interval formulas and fixed-weight effective units are correct", {
  f <- confidence_fixture()
  r <- confidence_make(f$model,f$samples)
  expected <- sqrt(log(2*8/.05)/(2*100))
  expect_equal(r$cells$lower,pmax(0,r$cells$estimate-expected))
  expect_equal(r$cells$upper,pmin(1,r$cells$estimate+expected))
  expect_equal(r$cells$effective_units,rep(100,8))
  expect_null(r$samples)
  for (n in c(1,2,10,100)) for (k in unique(c(0,1,n%/%2,n))) {
    z <- rep(c(0,1),c(n-k,k))
    r <- confidence_make(f$model,list(cell_1=z),"clopper_pearson")
    expect_equal(c(r$cells$lower[1],r$cells$upper[1]),
      as.numeric(stats::binom.test(k,n,conf.level=1-.05/8)$conf.int),tolerance=1e-12)
    expect_true(all(r$cells$lower[-1]==0&r$cells$upper[-1]==1))
  }
  r <- poc_confidence_region(f$model,list(cell_1=c(0,.5,1)),target="Fixed weighted cluster target",
    sampling_unit="cluster",weights=list(cell_1=c(1,2,3)),assume_sampling=TRUE)
  expect_equal(r$cells$estimate[1],2/3)
  expect_equal(r$cells$effective_units[1],36/14)
  expect_equal(r$cells$n_units[1],3)
  expect_equal(r$cells$lower[1],max(0,2/3-sqrt(14/36*log(16/.05)/2)))
  expect_identical(r$sampling_unit,"cluster")
  expect_equal(confidence_make(f$model,list(cell_1=numeric()))$cells$n_units,rep(0,8))
})

test_that("template numbers never become fixed confidence constraints", {
  f <- confidence_fixture()
  a <- matrix(c(.6,.4,.3,.7),2,byrow=TRUE,dimnames=dimnames(f$model$o))
  alternate <- poc_model(matrix(.25,2,2,dimnames=dimnames(a)),a)
  q <- poc_query(c("0"="0","1"="1"))
  first <- poc_confidence_bounds(confidence_make(f$model,f$samples),q)
  second <- poc_confidence_bounds(confidence_make(alternate,f$samples),q)
  expect_equal(first,second)
  # Empirical centers violate consistency, but a sufficiently broad region fits.
  s <- list(cell_1=rep(1,10),cell_5=rep(0,10))
  r <- confidence_make(f$model,s)
  expect_gt(r$cells$lower[1],r$cells$estimate[5])
  expect_identical(poc_confidence_bounds(r,q)$status,"ok")
  # Larger samples reject this same incompatible combination; do not drop it.
  s <- list(cell_1=rep(1,1000),cell_5=rep(0,1000))
  z <- poc_confidence_bounds(confidence_make(f$model,s),q)
  expect_identical(z$status,"empty_confidence_region")
  expect_true(is.na(z$lower)&&is.na(z$upper))
  failed_conditional <- poc_confidence_bounds(confidence_make(f$model,s),
    poc_query(c("0"="0"),"1","1",TRUE))
  expect_identical(failed_conditional$status,"empty_confidence_region")
  expect_true(is.na(failed_conditional$condition_may_be_zero))
  free <- poc_confidence_bounds(confidence_make(f$model,list()),q)
  expect_equal(c(free$lower,free$upper),c(0,1))
})

test_that("projection agrees with bands and independent endpoint calculations", {
  f <- confidence_fixture()
  r <- confidence_make(f$model,f$samples)
  delta <- sqrt(log(16/.05)/200)
  queries <- list(poc_query(c("0"="0","1"="1")),
    poc_query(c("0"="0"),observed_x="1",observed_y="1",conditional=TRUE),
    poc_query(c("1"="1"),observed_x="0",observed_y="0",conditional=TRUE),
    poc_query(c("0"="0","1"="1"),observed_y="1"))
  for (q in queries) {
    z <- poc_confidence_bounds(r,q,witnesses=TRUE)
    old <- poc_input_sensitivity(f$model,q,observed_delta=delta,interventional_delta=delta)
    expect_identical(z$status,"ok")
    expect_equal(c(z$lower,z$upper),c(old$lower,old$upper),tolerance=1e-8)
    for (side in c("lower","upper")) {
      w <- attr(z,"witnesses")[[side]]
      mass <- w$probability
      expect_equal(sum(mass),1,tolerance=1e-8)
      expect_gte(min(mass),-1e-8)
      ox <- w$factual
      oy <- w$response[cbind(seq_along(mass),ox)]
      margins <- c(unlist(lapply(1:2,function(x) vapply(1:2,function(y)
        sum(mass[ox==x&oy==y]),0))),unlist(lapply(1:2,function(x)
        vapply(1:2,function(y) sum(mass[w$response[,x]==y]),0))))
      expect_true(all(margins>=r$cells$lower-1e-8&margins<=r$cells$upper+1e-8))
      event <- rep(TRUE,length(mass))
      for (x in names(q$counterfactual)) event <- event &
        w$response[,match(x,rownames(f$model$o))]==match(q$counterfactual[[x]],colnames(f$model$o))
      cond <- rep(TRUE,length(mass))
      if(!is.null(q$observed_x)) cond <- cond & ox==match(q$observed_x,rownames(f$model$o))
      if(!is.null(q$observed_y)) cond <- cond & oy==match(q$observed_y,colnames(f$model$o))
      truth <- sum(mass[event&cond]) / if(q$conditional) sum(mass[cond]) else 1
      expect_equal(truth,z[[side]],tolerance=1e-8)
    }
  }
  # With no intervention inputs, binary PNS maximum has an independent formula.
  r <- confidence_make(f$model,f$samples[1:4],method="clopper_pearson")
  z <- poc_confidence_bounds(r,queries[[1]])
  upper <- min(1,r$cells$upper[1]+r$cells$upper[4],1-r$cells$lower[2]-r$cells$lower[3])
  expect_equal(c(z$lower,z$upper),c(0,upper),tolerance=1e-8)
  # A deterministic large-sample fixture checks nonvacuous conditional bounds.
  # Repeated values are a software fixture, not evidence of independent papers.
  precise <- confidence_make(f$model,lapply(f$samples,rep,times=100))
  q <- queries[[2]]
  z <- poc_confidence_bounds(precise,q)
  identified <- poc_exact(f$model,q)
  expect_gt(z$lower,.5)
  expect_lte(z$lower,identified$lower+1e-8)
  expect_gte(z$upper,identified$upper-1e-8)
  expect_false(z$condition_may_be_zero)
  expect_gt(z$condition_probability_min,0)
})

test_that("one-period static and longitudinal confidence projections agree", {
  f <- confidence_fixture()
  o <- expand.grid(history=c("0","1"),path=c("0","1"),stringsAsFactors=FALSE)
  o$prob <- as.vector(f$model$o)
  a <- data.frame(regime=o$history,path=o$path,prob=as.vector(f$model$a))
  # Map intervention rows in row-major order to match the static map.
  a <- a[c(1,3,2,4),]
  model <- poc_longitudinal_model(o,matrix(0:1,ncol=1,dimnames=list(c("0","1"),"t1")),a)
  for (method in c("hoeffding","clopper_pearson")) {
    r <- confidence_make(model,f$samples,method)
    sr <- confidence_make(f$model,f$samples,method)
    sq <- list(pns=poc_query(c("0"="0","1"="1")),
      pn=poc_query(c("0"="0"),"1","1",TRUE),ps=poc_query(c("1"="1"),"0","0",TRUE))
    for (kind in names(sq)) {
      z <- poc_confidence_bounds(r,poc_longitudinal_query(kind,"1","0",1))
      s <- poc_confidence_bounds(sr,sq[[kind]])
      expect_identical(z$status,s$status)
      expect_equal(c(z$lower,z$upper,z$condition_probability_min,z$condition_probability_max),
        c(s$lower,s$upper,s$condition_probability_min,s$condition_probability_max),tolerance=1e-8)
    }
  }
})

test_that("zero conditions and confidence labels survive reporting", {
  f <- confidence_fixture()
  r <- confidence_make(f$model,list())
  q <- poc_query(c("0"="0"),"1","1",TRUE)
  z <- poc_confidence_bounds(r,q)
  expect_true(z$condition_may_be_zero)
  expect_equal(c(z$condition_probability_min,z$condition_probability_max),c(0,1))
  expect_true(is.na(z$denominator))
  expect_equal(c(z$lower,z$upper),c(0,1))
  # Outcome state is structurally excluded by allowed_paths, not a sample zero.
  o <- data.frame(history=c("a","b"),path="0",prob=.5)
  model <- poc_longitudinal_model(o,matrix(0:1,ncol=1,dimnames=list(c("a","b"),"t")),allowed_paths="0")
  undefined <- poc_confidence_bounds(confidence_make(model,list()),poc_longitudinal_query("pn","b","a",1))
  expect_identical(undefined$status,"undefined_condition")
  expect_equal(undefined$condition_probability_max,0)
  expect_identical(summary(z)$bounds$uncertainty_kind,"simultaneous_identified_set_confidence_region")
  path <- tempfile()
  poc_export_summary(summary(z),path)
  saved <- read.delim(file.path(path,"bounds.tsv"),check.names=FALSE)
  expect_equal(saved$conf_level,.95)
  expect_identical(saved$target,r$target)
  expect_identical(saved$region_method,"hoeffding_bonferroni")
  unlink(path,recursive=TRUE)
})

test_that("multivalued projections preserve event families and static reductions", {
  types <- expand.grid(y1=1:3,y2=1:3,y3=1:3,x=1:3)
  mass <- seq_len(nrow(types)); N <- sum(mass)
  o <- a <- matrix(0,3,3,dimnames=list(c("g1","g2","g3"),c("low","middle","high")))
  response <- as.matrix(types[,1:3])
  actual <- response[cbind(seq_len(nrow(types)),types$x)]
  for(x in 1:3) for(y in 1:3) {
    o[x,y] <- sum(mass[types$x==x&actual==y])/N
    a[x,y] <- sum(mass[response[,x]==y])/N
  }
  model <- poc_model(o,a)
  samples <- setNames(lapply(c(as.vector(t(o)),as.vector(t(a))),function(p)
    rep(c(0,1),c(N-round(N*p),round(N*p)))),poc_margin_map(model)$cell_id)
  r <- confidence_make(model,samples)
  delta <- sqrt(log(36/.05)/(2*N))
  general <- poc_query(c(g1="low",g2="middle",g3="high"),observed_y="middle",conditional=TRUE)
  z <- poc_confidence_bounds(r,general)
  old <- poc_input_sensitivity(model,general,observed_delta=delta,interventional_delta=delta)
  expect_identical(z$status,"ok")
  expect_equal(c(z$lower,z$upper),c(old$lower,old$upper),tolerance=1e-8)
  rows <- expand.grid(outcome=colnames(o),treatment=rownames(o),stringsAsFactors=FALSE)
  observed <- data.frame(history=rows$treatment,path=rows$outcome,prob=as.vector(t(o)))
  margins <- data.frame(regime=rows$treatment,path=rows$outcome,prob=as.vector(t(a)))
  paths <- matrix(colnames(o),ncol=1,dimnames=list(colnames(o),"t1"))
  regimes <- matrix(1:3,ncol=1,dimnames=list(rownames(o),"t1"))
  longitudinal <- poc_longitudinal_model(observed,regimes,margins,outcome_paths=paths)
  lr <- confidence_make(longitudinal,samples)
  for(active in colnames(o)) for(reference in colnames(o)) for(kind in c("pns","pn","ps")) {
    sq <- switch(kind,pns=poc_query(c(g3=active,g1=reference)),
      pn=poc_query(c(g1=reference),"g3",active,TRUE),ps=poc_query(c(g3=active),"g1",reference,TRUE))
    lq <- poc_longitudinal_query(kind,"g3","g1",1,regime_outcome=active,reference_outcome=reference)
    s <- poc_confidence_bounds(r,sq); l <- poc_confidence_bounds(lr,lq)
    expect_identical(s$status,l$status)
    expect_equal(c(s$lower,s$upper),c(l$lower,l$upper),tolerance=1e-8)
  }
})

test_that("unsafe sampling inputs and unsupported templates are rejected", {
  f <- confidence_fixture()
  expect_error(poc_margin_map(list()),"template")
  expect_error(poc_confidence_region(f$model,list(),target="x",sampling_unit="individual"),"assume_sampling")
  expect_error(poc_confidence_region(f$model,list(),target="x",assume_sampling=TRUE),"sampling_unit")
  expect_error(confidence_make(f$model,list(unknown=1)),"identifiers")
  expect_error(confidence_make(f$model,list(1)),"named list")
  expect_error(confidence_make(f$model,list(cell_1=1,cell_1=0)),"named list")
  for (z in list(NA_real_,Inf,-1,2,"0",matrix(0),list(0)))
    expect_error(confidence_make(f$model,list(cell_1=z)),"numeric vector")
  expect_error(confidence_make(f$model,list(cell_1=.5),"clopper_pearson"),"binary")
  expect_error(confidence_make(f$model,list(cell_1=0),"clopper_pearson",weights=list()),"unweighted")
  expect_error(poc_confidence_region(f$model,list(),target="x",sampling_unit="cluster",
    method="clopper_pearson",assume_sampling=TRUE),"unweighted")
  for (w in list(c(0,0),c(-1,2),NA_real_,c(1,Inf),1,character()))
    expect_error(confidence_make(f$model,list(cell_1=c(0,1)),weights=list(cell_1=w)),"Weights")
  expect_error(confidence_make(f$model,list(),weights=list(cell_1=1)),"supplied sample")
  for (level in c(0,1,-1,Inf,NA_real_))
    expect_error(confidence_make(f$model,list(),conf_level=level),"confidence level")
  r <- confidence_make(f$model,list())
  expect_error(poc_confidence_bounds(r,poc_query(c("0"="0")),max_variables=1),"limit")
  expect_error(poc_confidence_bounds(r,poc_query(c("0"="0")),witnesses=NA),"flag")
  r$cells <- r$cells[nrow(r$cells):1,]
  expect_error(poc_confidence_bounds(r,poc_query(c("0"="0"))),"no longer match")
})
