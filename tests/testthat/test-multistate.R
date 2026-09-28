multistate_fixture <- function(periods = 2L, groups = 2L) {
  levels <- c("low", "middle", "high")
  paths <- as.matrix(expand.grid(rep(list(levels), periods), stringsAsFactors = FALSE))
  rownames(paths) <- paste0("path_", seq_len(nrow(paths)))
  regimes <- matrix(rep(seq_len(groups)-1L, each = periods), groups, periods, byrow = TRUE,
                    dimnames = list(paste0("g", seq_len(groups)), paste0("t",seq_len(periods))))
  types <- as.matrix(expand.grid(rep(list(seq_len(nrow(paths))), groups)))
  response <- types[rep(seq_len(nrow(types)), groups), , drop = FALSE]
  factual <- rep(seq_len(groups), each = nrow(types))
  mass <- (seq_len(nrow(response)) %% 11)+1
  mass <- mass/sum(mass)
  actual <- response[cbind(seq_len(nrow(response)), factual)]
  obs <- expand.grid(history = rownames(regimes), path = rownames(paths), stringsAsFactors = FALSE)
  obs$prob <- vapply(seq_len(nrow(obs)), function(i)
    sum(mass[factual == match(obs$history[i], rownames(regimes)) &
               actual == match(obs$path[i], rownames(paths))]), 0.0)
  margins <- obs
  names(margins)[1L] <- "regime"
  margins$prob <- vapply(seq_len(nrow(margins)), function(i)
    sum(mass[response[, match(margins$regime[i], rownames(regimes))] ==
               match(margins$path[i], rownames(paths))]), 0.0)
  times <- expand.grid(regime = rownames(regimes), horizon = seq_len(periods),
                       outcome = levels, stringsAsFactors = FALSE)
  times$prob <- vapply(seq_len(nrow(times)), function(i)
    sum(mass[paths[response[, match(times$regime[i], rownames(regimes))], times$horizon[i]] ==
               times$outcome[i]]), 0.0)
  list(paths=paths,regimes=regimes,observed=obs,margins=margins,times=times,
       response=response,factual=factual,mass=mass)
}

test_that("multivalued one-period queries reduce to the static exact LP", {
  f <- multistate_fixture(1L,3L)
  model <- poc_longitudinal_model(f$observed,f$regimes,time_margins=f$times,
                                  outcome_paths=f$paths)
  paths_model <- poc_longitudinal_model(f$observed,f$regimes,f$margins,
                                        outcome_paths=f$paths)
  o <- a <- matrix(0,3L,3L,dimnames=list(rownames(f$regimes),as.character(f$paths[,1L])))
  for (i in seq_len(nrow(f$observed)))
    o[f$observed$history[i], f$paths[f$observed$path[i],1L]] <- f$observed$prob[i]
  for (i in seq_len(nrow(f$times))) a[f$times$regime[i],f$times$outcome[i]] <- f$times$prob[i]
  static <- poc_model(o,a)
  for (g in rownames(o)) for (r in setdiff(rownames(o),g))
    for (yg in colnames(o)) for (yr in colnames(o)) for (kind in c("pns","pn","ps")) {
      q <- poc_longitudinal_query(kind,g,r,1,regime_outcome=yg,reference_outcome=yr)
      sq <- switch(kind,
        pns=poc_query(setNames(c(yg,yr),c(g,r))),
        pn=poc_query(setNames(yr,r),observed_x=g,observed_y=yg,conditional=TRUE),
        ps=poc_query(setNames(yg,g),observed_x=r,observed_y=yr,conditional=TRUE))
      expected <- poc_exact(static,sq)
      for (m in list(model,paths_model)) {
        actual <- poc_longitudinal_bounds(m,q)
        expect_equal(actual$lower,expected$lower,tolerance=1e-9)
        expect_equal(actual$upper,expected$upper,tolerance=1e-9)
        expect_equal(actual$denominator,expected$denominator,tolerance=1e-9)
        expect_identical(actual$status,"ok")
      }
    }
})

test_that("categorical path events and unions contain directly computed truth", {
  f <- multistate_fixture()
  model <- poc_longitudinal_model(f$observed,f$regimes,f$margins,outcome_paths=f$paths)
  active <- c("middle","high"); reference <- "low"
  a <- matrix(f$paths[f$response[,2L],] %in% active,nrow(f$response),2L)
  b <- matrix(f$paths[f$response[,1L],] %in% reference,nrow(f$response),2L)
  for (h in 1:2) for (kind in c("pns","pn","ps","persistent","event_time","lagged")) {
    if (kind=="lagged" && h==1L) next
    q <- poc_longitudinal_query(kind,"g2","g1",h,regime_outcome=active,reference_outcome=reference)
    num <- switch(kind,pns=a[,h]&b[,h],pn=b[,h],ps=a[,h],
      persistent=rowSums(a[,seq_len(h),drop=FALSE]&b[,seq_len(h),drop=FALSE])==h,
      event_time=a[,h]&b[,h]&(if(h==1L) TRUE else !a[,1L]),
      lagged=a[,2L]&b[,1L])
    given <- switch(kind,pn=f$factual==2L&a[,h],ps=f$factual==1L&b[,h],rep(TRUE,length(num)))
    truth <- sum(f$mass[num&given])/sum(f$mass[given])
    result <- poc_longitudinal_bounds(model,q)
    expect_lte(result$lower,truth+1e-9)
    expect_gte(result$upper,truth-1e-9)
    expect_equal(result$denominator,sum(f$mass[given]),tolerance=1e-9)
    expect_identical(result$regime_outcome,"middle;high")
    expect_identical(result$reference_outcome,"low")
    zero <- poc_longitudinal_sensitivity(model,q)
    wide <- poc_longitudinal_sensitivity(model,q,observed_delta=.05,path_delta=.1)
    expect_equal(c(zero$lower,zero$upper),c(result$lower,result$upper),tolerance=1e-9)
    expect_lte(wide$lower,truth+1e-9)
    expect_gte(wide$upper,truth-1e-9)
  }
  timed <- poc_longitudinal_model(f$observed,f$regimes,time_margins=f$times,outcome_paths=f$paths)
  q <- poc_longitudinal_query("pns","g2","g1",2,regime_outcome=active,reference_outcome=reference)
  exact <- poc_longitudinal_bounds(timed,q)
  zero <- poc_longitudinal_sensitivity(timed,q)
  wide <- poc_longitudinal_sensitivity(timed,q,time_delta=.1)
  expect_equal(c(zero$lower,zero$upper),c(exact$lower,exact$upper),tolerance=1e-9)
  expect_lte(wide$lower,exact$lower+1e-9)
  expect_gte(wide$upper,exact$upper-1e-9)
})

test_that("all categorical point-mass path pairs identify all seven query families", {
  f <- multistate_fixture()
  for (p0 in rownames(f$paths)) for (p1 in rownames(f$paths)) {
    observed <- f$observed
    observed$prob <- ifelse(observed$history=="g1" & observed$path==p0,.4,
      ifelse(observed$history=="g2" & observed$path==p1,.6,0))
    margins <- f$margins
    margins$prob <- as.numeric((margins$regime=="g1" & margins$path==p0) |
                                (margins$regime=="g2" & margins$path==p1))
    model <- poc_longitudinal_model(observed,f$regimes,margins,outcome_paths=f$paths)
    a <- f$paths[p1,] %in% c("middle","high")
    b <- f$paths[p0,]=="low"
    for (kind in c("pns","pn","ps","persistent","event_time","lagged","trajectory")) {
      q <- if(kind=="trajectory") poc_longitudinal_query(kind,"g2","g1",2,
        regime_path=p1,reference_path=p0) else poc_longitudinal_query(kind,"g2","g1",2,
        regime_outcome=c("middle","high"),reference_outcome="low")
      result <- poc_longitudinal_bounds(model,q)
      undefined <- (kind=="pn" && !a[2L]) || (kind=="ps" && !b[2L])
      if (undefined) {
        expect_identical(result$status,"undefined_condition")
        expect_true(is.na(result$lower)&&is.na(result$upper))
      } else {
        truth <- switch(kind,pns=a[2L]&b[2L],pn=b[2L],ps=a[2L],
          persistent=all(a&b),event_time=!a[1L]&a[2L]&b[2L],lagged=a[2L]&b[1L],trajectory=TRUE)
        expect_identical(result$status,"ok")
        expect_equal(c(result$lower,result$upper),rep(as.numeric(truth),2L),tolerance=1e-9)
      }
    }
  }
})

test_that("ordered paths, full categorical no anticipation and prefix sets are enforced", {
  f <- multistate_fixture()
  levels <- c("low","middle","high")
  keep <- match(f$paths[,1L],levels)<=match(f$paths[,2L],levels)
  p <- f$paths[keep,,drop=FALSE]
  o <- expand.grid(history=c("off","late","early"),path=rownames(p),stringsAsFactors=FALSE)
  regimes <- rbind(off=c(0,0),late=c(0,1),early=c(1,1))
  o$prob <- ifelse(o$path==rownames(p)[1L],1/3,0)
  model <- poc_longitudinal_model(o,regimes,outcome_paths=p,outcome_order=levels,absorbing=TRUE)
  expect_identical(model$outcome_order,levels)
  # At t1 all are low. PN conditioned on {low,middle} pools late+off.
  q <- poc_longitudinal_query("pn","late","early",1,condition_on="prefix",
    regime_outcome=c("low","middle"),reference_outcome="low")
  z <- poc_longitudinal_bounds(model,q)
  expect_equal(z$denominator,2/3)
  expect_identical(z$outcome_order,"low;middle;high")
  expect_equal(poc_longitudinal_sensitivity(model,q)$denominator,2/3)
  # Shared prefix cannot be low under off and middle under late for everyone.
  margins <- data.frame(regime=c("off","late"),horizon=1,
    outcome=c("low","middle"),prob=1)
  conflicting_obs <- o
  conflicting_obs$prob[conflicting_obs$history=="late"] <- 0
  middle_path <- rownames(p)[p[,1L]=="middle" & p[,2L]=="middle"]
  conflicting_obs$prob[conflicting_obs$history=="late" & conflicting_obs$path==middle_path] <- 1/3
  expect_error(poc_longitudinal_model(conflicting_obs,regimes,time_margins=margins,outcome_paths=p),"No finite")
  relaxed <- poc_longitudinal_model(conflicting_obs,regimes,time_margins=margins,
    outcome_paths=p,no_anticipation=FALSE)
  expect_s3_class(relaxed,"pnsbib_longitudinal_model")
})

test_that("multistate malformed domains, orders, margins and queries are rejected", {
  f <- multistate_fixture()
  make <- function(p=f$paths,...) poc_longitudinal_model(f$observed,f$regimes,outcome_paths=p,...)
  expect_error(make(unname(f$paths)),"outcome_paths")
  expect_error(make(f$paths[,1L,drop=FALSE]),"outcome_paths")
  bad <- f$paths; bad[1,1] <- NA
  expect_error(make(bad),"outcome_paths")
  bad <- f$paths; bad[1,] <- bad[2,]
  expect_error(make(bad),"distinct")
  expect_error(make(absorbing=TRUE),"outcome_order")
  expect_error(make(absorbing=TRUE,outcome_order=c("low","middle","high")),"violate")
  expect_error(make(outcome_order=c("low","high")),"every path state")
  expect_error(make(outcome_order=c("low","low","middle","high")),"every path state")
  expect_error(make(max_variables=10),"max_variables")
  bad <- f$times; bad$prob1 <- bad$prob
  expect_error(make(time_margins=bad),"schemas")
  expect_error(make(time_margins=f$times[,c("regime","horizon","prob")]),"schemas")
  expect_error(make(time_margins=rbind(f$times,f$times[1L,])),"time margins")
  bad <- f$times; bad$outcome[1L] <- "absent"
  expect_error(make(time_margins=bad),"time margins")
  bad <- f$times; bad$prob[1L] <- Inf
  expect_error(make(time_margins=bad),"time margins")
  for (invalid_time in c(NA_real_, NaN, Inf, -Inf, 1e20, .5)) {
    bad <- f$times; bad$horizon[1L] <- invalid_time
    expect_error(make(time_margins=bad),"time margins")
    expect_error(poc_longitudinal_query("pns","g2","g1",invalid_time),"time indices")
    expect_error(poc_longitudinal_query("pns","g2","g1",2,start=invalid_time),"time indices")
    expect_error(poc_longitudinal_query("lagged","g2","g1",2,lag=invalid_time),"time indices")
  }
  model <- make()
  expect_error(poc_longitudinal_bounds(model,poc_longitudinal_query("pns","g2","g1",1)),"outcome set")
  expect_error(poc_longitudinal_query("pns","g2","g1",1,regime_outcome=c("high","high")),"Outcome sets")
  expect_error(poc_longitudinal_query("trajectory","g2","g1",1,regime_path="p1",
    reference_path="p2",regime_outcome="high"),"inapplicable")
  # Declared but structurally excluded states remain legal zero events.
  reduced <- poc_longitudinal_model(f$observed[f$observed$path=="path_1",] |>
    transform(prob=c(.4,.6)),f$regimes,outcome_paths=f$paths,allowed_paths="path_1")
  expect_identical(poc_longitudinal_bounds(reduced,poc_longitudinal_query("pn","g2","g1",1,
    regime_outcome="high",reference_outcome="low"))$status,"undefined_condition")
})
