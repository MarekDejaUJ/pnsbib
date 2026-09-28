curve_fixture <- function() {
  regimes <- rbind(off=c(0,0),on=c(1,1))
  observed <- expand.grid(history=c("off","on"),path=c("00","01","11"),stringsAsFactors=FALSE)
  observed$prob <- c(.2,.1,.1,.1,.2,.3)
  poc_longitudinal_model(observed,regimes,absorbing=TRUE)
}

test_that("curves preserve explicit grid order and exact single-call results", {
  model <- curve_fixture()
  result <- poc_curve(model,"on","off",horizons=c(2,1),kind=c("ps","pns","pn"),condition_on="prefix")
  expect_identical(result$horizon,c(2L,2L,2L,1L,1L,1L))
  expect_identical(result$estimand,rep(c("PS","PNS","PN"),2))
  expect_identical(result$curve_row,1:6)
  for (i in 1:6) {
    kind <- tolower(result$estimand[i])
    q <- poc_longitudinal_query(kind,"on","off",result$horizon[i],
      condition_on=if(kind=="pns")"history" else "prefix")
    single <- poc_longitudinal_bounds(model,q)
    row <- result[i,names(single),drop=FALSE]; rownames(row) <- NULL
    expect_equal(row,single)
  }
  expect_s3_class(summary(result),"summary.pnsbib_analysis")
})

test_that("curve bands, floors, categorical labels and exports preserve their meaning", {
  model <- curve_fixture()
  exact <- poc_curve(model,"on","off")
  zero <- poc_curve(model,"on","off",bands=list())
  expect_equal(zero$lower,exact$lower,tolerance=1e-9)
  expect_equal(zero$upper,exact$upper,tolerance=1e-9)
  wide <- poc_curve(model,"on","off",bands=list(observed_delta=.05,
    minimum_condition_probability=.1))
  expect_true(all(is.na(wide$minimum_condition_probability[wide$estimand=="PNS"])))
  expect_equal(wide$minimum_condition_probability[wide$estimand!="PNS"],rep(.1,4))
  expect_true(all(wide$lower<=exact$lower+1e-9))
  expect_true(all(wide$upper>=exact$upper-1e-9))
  paths <- matrix(c("low","low","low","high","high","high"),3,2,byrow=TRUE,
                  dimnames=list(c("00","01","11"),NULL))
  obs <- expand.grid(history=c("off","on"),path=rownames(paths),stringsAsFactors=FALSE)
  obs$prob <- c(.2,.1,.1,.1,.2,.3)
  multistate <- poc_longitudinal_model(obs,model$regimes,outcome_paths=paths,
    absorbing=TRUE,outcome_order=c("low","high"))
  same <- poc_curve(multistate,"on","off",regime_outcome="high",reference_outcome="low")
  expect_equal(same$lower,exact$lower,tolerance=1e-9)
  expect_equal(same$upper,exact$upper,tolerance=1e-9)
  expect_identical(same$regime_outcome,rep("high",6))
  out <- tempfile("curve-export-")
  poc_export_summary(poc_analysis(same,provenance=list(test="categorical_curve")),out)
  back <- read.delim(file.path(out,"bounds.tsv"))
  expect_equal(back$lower,same$lower,tolerance=1e-12)
  expect_identical(back$outcome_order,same$outcome_order)
  expect_identical(back$regime_outcome,same$regime_outcome)
})

test_that("curve statuses and nonmonotone PNS are not altered", {
  model <- curve_fixture()
  # Active always crosses at t1, reference catches up at t2.
  a <- data.frame(regime=c("off","on"),path=c("01","11"),prob=1)
  obs <- expand.grid(history=c("off","on"),path=c("00","01","11"),stringsAsFactors=FALSE)
  obs$prob <- c(0,0,.5,0,0,.5)
  model <- poc_longitudinal_model(obs,model$regimes,a,absorbing=TRUE)
  curve <- poc_curve(model,"on","off")
  expect_equal(curve$lower[curve$estimand=="PNS"],c(1,0))
  expect_equal(curve$upper[curve$estimand=="PNS"],c(1,0))
  expect_identical(curve$status[curve$estimand=="PS"],c("ok","undefined_condition"))
  # Serialized binary models/queries from <=0.0.5 still work.
  legacy <- model; legacy$outcome_order <- legacy$outcome_levels <- NULL
  query <- poc_longitudinal_query("pns","on","off",1)
  query$regime_outcome <- query$reference_outcome <- NULL
  expect_equal(poc_longitudinal_bounds(legacy,query)$lower,1)
  expect_equal(poc_longitudinal_sensitivity(legacy,query)$lower,1)
})

test_that("invalid curve grids and sensitivity options fail explicitly", {
  m <- curve_fixture()
  expect_error(poc_curve(list(),"on","off"),"model")
  for (h in list(numeric(),c(1,1),0,3,NA,Inf,1.5))
    expect_error(poc_curve(m,"on","off",horizons=h),"horizons")
  for (k in list(character(),c("pn","pn"),"trajectory","typo"))
    expect_error(poc_curve(m,"on","off",kind=k),"kinds")
  expect_error(poc_curve(m,"on","off",kind="lagged"),"lag")
  expect_error(poc_curve(m,"on","off",bands=list(.1)),"named list")
  expect_error(poc_curve(m,"on","off",bands=list(unknown=.1)),"named list")
  expect_error(poc_curve(m,"on","off",bands=list(time_delta=-.1)),"time_delta")
  expect_error(poc_curve(m,"on","off",kind="pns",bands=list(minimum_condition_probability=.1)),"floor")
  expect_s3_class(poc_curve(m,"on","off",horizons=2,kind=c("persistent","event_time","lagged")),
    "pnsbib_longitudinal_result")
})
