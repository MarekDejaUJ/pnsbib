shu_fixture <- function(m = 3L, n = 3L, weights = NULL) {
  grid <- as.matrix(expand.grid(c(list(x = seq_len(m)),
    setNames(rep(list(seq_len(n)), m), paste0("r", seq_len(m))))))
  if (is.null(weights)) weights <- seq_len(nrow(grid)) %% 11 + 1
  weights <- weights / sum(weights)
  o <- a <- matrix(0, m, n, dimnames = list(paste0("x", seq_len(m)), paste0("y", seq_len(n))))
  factual <- grid[cbind(seq_len(nrow(grid)), grid[, 1L] + 1L)]
  for (x in seq_len(m)) for (y in seq_len(n)) {
    o[x,y] <- sum(weights[grid[,1L] == x & factual == y])
    a[x,y] <- sum(weights[grid[,x+1L] == y])
  }
  list(model = poc_model(o,a), grid = grid, weights = weights)
}

test_that("Shu closed forms reproduce both source examples and corrected singletons", {
  tables <- list(
    list(o=rbind(c(131,68,1),c(45,22,51),c(38,483,61))/900,
         a=rbind(c(46,23,231),c(270,8,22),c(40,223,37))/300,
         q=poc_query(c(x1="y3",x2="y1",x3="y2")), interval=c(458/900,529/900)),
    list(o=rbind(c(67,129,193),c(11,17,87),c(53,53,70),c(46,436,38))/1200,
         a=rbind(c(195,51,54),c(11,266,23),c(80,198,22),c(100,147,53))/300,
         q=poc_query(c(x1="y1",x2="y2",x3="y2"),"x4","y2",TRUE), interval=c(15/436,1)))
  for (z in tables) {
    dimnames(z$o) <- dimnames(z$a) <- list(paste0("x",seq_len(nrow(z$o))),paste0("y",1:3))
    fit <- poc_shu2026(poc_model(z$o,z$a),z$q)
    expect_equal(c(fit$lower,fit$upper),z$interval,tolerance=1e-10)
    expect_identical(fit$sharpness,"valid_closed_form")
  }
  o <- rbind(c(.1,.1,.3),c(.1,.125,.025),c(.1,.125,.025))
  a <- rbind(c(.55,.125,.325),c(.35,.375,.275),c(.35,.375,.275))
  dimnames(o) <- dimnames(a) <- list(paste0("x",1:3),paste0("y",1:3))
  model <- poc_model(o,a)
  for (y in colnames(o)) {
    q <- poc_query(c(x1="y1"),observed_y=y)
    closed <- poc_shu2026(model,q); exact <- poc_exact(model,q)
    expect_equal(c(closed$lower,closed$upper),c(exact$lower,exact$upper),tolerance=1e-9)
  }
})

test_that("distinct subset optimization agrees with explicit enumeration", {
  for (d in list(c(.01,.2),c(.7,.1,.4),c(0,.3,.3,.9),seq(0,1,length.out=9))) {
    expected <- vapply(2:length(d),function(k)
      min(combn(d,k,FUN=sum))/(k-1),0.0)
    expect_equal(pnsbib:::.poc_shu_subsets(d),expected)
    expect_equal(pnsbib:::.poc_shu_subsets(rev(d)),expected)
  }
  expect_length(pnsbib:::.poc_shu_subsets(numeric()),0)
  expect_length(pnsbib:::.poc_shu_subsets(.1),0)
})

test_that("repeated factual matches retain their disjoint-treatment subset bound", {
  o <- rbind(x1=c(y1=.4,y2=.1),x2=c(y1=.4,y2=.1))
  a <- rbind(x1=c(y1=.5,y2=.5),x2=c(y1=.5,y2=.5))
  model <- poc_model(o,a)
  query <- poc_query(c(x1="y1",x2="y1"),observed_y="y1")
  result <- poc_shu2026(model,query)
  exact <- poc_exact(model,query)
  # Using just C=.8 and a_i=.5 would miss the subset upper bound .1+.1=.2.
  expect_equal(c(result$lower,result$upper),c(0,.2),tolerance=1e-9)
  expect_equal(c(result$lower,result$upper),c(exact$lower,exact$upper),tolerance=1e-9)
  expect_match(result$method,"project repeated-outcome extension",fixed=TRUE)
  query$conditional <- TRUE
  expect_equal(poc_shu2026(model,query)$upper,.25,tolerance=1e-9)
})

test_that("every binary point mass identifies every closed-form event", {
  base <- shu_fixture(2,2)
  qs <- as.matrix(expand.grid(rep(list(0:2),2)))[-1,,drop=FALSE]
  for (row in seq_len(nrow(base$grid))) {
    z <- shu_fixture(2,2,as.numeric(seq_len(nrow(base$grid))==row))
    for (v in seq_len(nrow(qs))) for (ox in 0:2) for (oy in 0:2) {
      x <- which(qs[v,]>0); y <- qs[v,x]
      q <- poc_query(setNames(paste0("y",y),paste0("x",x)),
        if(ox) paste0("x",ox) else NULL,if(oy) paste0("y",oy) else NULL)
      truth <- as.numeric(all(z$grid[row,x+1L]==y) &&
        (!ox || z$grid[row,1L]==ox) &&
        (!oy || z$grid[row,z$grid[row,1L]+1L]==oy))
      fit <- poc_shu2026(z$model,q)
      expect_equal(c(fit$lower,fit$upper),rep(truth,2),tolerance=1e-9)
      if (ox || oy) {
        q$conditional <- TRUE
        cond <- poc_shu2026(z$model,q)
        den <- if (ox && oy) z$model$o[ox,oy] else if(ox) z$model$px[ox] else z$model$py[oy]
        expect_equal(cond$denominator,unname(den))
        if (den == 0) expect_identical(cond$status,"undefined_condition") else
          expect_equal(c(cond$lower,cond$upper),rep(truth,2))
      }
    }
  }
})

test_that("four families and repeated labels contain the same-information exact interval", {
  for (dims in list(c(2,2),c(3,3),c(4,3),c(3,4))) {
    model <- shu_fixture(dims[1],dims[2])$model
    vectors <- list(c(x1="y1"),c(x1="y1",x2="y1"),c(x1="y1",x2="y2"),
      setNames(rep("y1",dims[1]),rownames(model$o)))
    for (cf in vectors) for (ox in c(NA_integer_,seq_len(dims[1])))
      for (oy in c(NA_integer_,seq_len(dims[2]))) {
        query <- poc_query(cf,if(is.na(ox)) NULL else rownames(model$o)[ox],
                           if(is.na(oy)) NULL else colnames(model$o)[oy])
        closed <- poc_shu2026(model,query); exact <- poc_exact(model,query)
        expect_lte(closed$lower,exact$lower+1e-8)
        expect_gte(closed$upper,exact$upper-1e-8)
        expect_identical(closed$sharpness,"valid_closed_form")
        query$counterfactual <- rev(cf)
        reverse <- poc_shu2026(model,query)
        expect_equal(c(reverse$lower,reverse$upper),c(closed$lower,closed$upper))
        if (!is.na(ox) || !is.na(oy)) {
          query$conditional <- TRUE
          conditional <- poc_shu2026(model,query)
          expect_equal(c(conditional$lower,conditional$upper),
            c(closed$lower,closed$upper)/conditional$denominator)
        }
      }
    q <- poc_query(setNames(rep("y1",dims[1]),rownames(model$o)))
    joint <- poc_shu2026(model,q)
    q$observed_y <- "y1"
    redundant <- poc_shu2026(model,q)
    expect_match(redundant$method,"project repeated-outcome extension",fixed=TRUE)
    expect_equal(c(joint$lower,joint$upper),c(redundant$lower,redundant$upper))
  }
})

test_that("rectangular padding and relabelling preserve the event and bounds", {
  for (dims in list(c(2,3),c(3,2))) {
    model <- shu_fixture(dims[1],dims[2])$model
    if (dims[1] < dims[2]) {
      o <- rbind(model$o,pad=0)
      a <- rbind(model$a,pad=c(1,rep(0,dims[2]-1)))
    } else {
      o <- cbind(model$o,pad=0); a <- cbind(model$a,pad=0)
    }
    padded <- poc_model(o,a)
    for (cf in list(c(x1="y1"),c(x1="y1",x2="y1"),c(x1="y1",x2="y2")))
      for (ox in c(NA_character_,rownames(model$o))) for(oy in c(NA_character_,colnames(model$o))) {
        query <- poc_query(cf,if(is.na(ox)) NULL else ox,if(is.na(oy)) NULL else oy)
        before <- poc_shu2026(model,query); after <- poc_shu2026(padded,query)
        e1 <- poc_exact(model,query); e2 <- poc_exact(padded,query)
        expect_equal(c(before$lower,before$upper),c(after$lower,after$upper))
        expect_equal(c(e1$lower,e1$upper),c(e2$lower,e2$upper),tolerance=1e-8)
      }
    # Permute axes and rename them independently; arbitrary labels have no order.
    rx <- rev(seq_len(nrow(model$o))); ry <- rev(seq_len(ncol(model$o)))
    o <- model$o[rx,ry]; a <- model$a[rx,ry]
    dimnames(o) <- dimnames(a) <- list(paste0("policy-",rx),paste0("state-",ry))
    renamed <- poc_model(o,a)
    before <- poc_shu2026(model,poc_query(c(x1="y2",x2="y1"),observed_y="y1"))
    after <- poc_shu2026(renamed,poc_query(c("policy-1"="state-2","policy-2"="state-1"),observed_y="state-1"))
    expect_equal(c(before$lower,before$upper),c(after$lower,after$upper))
  }
})

test_that("closed-form API preserves denominator, status and method in exports", {
  model <- shu_fixture()$model
  q <- poc_query(c(x1="y1",x2="y2"),observed_x="x1",conditional=TRUE)
  result <- poc_shu2026(model,q); exact <- poc_exact(model,q)
  expect_equal(result$denominator,unname(model$px[1]))
  expect_equal(c(result$lower,result$upper),c(exact$lower,exact$upper),tolerance=1e-8)
  q$observed_y <- "y3"
  expect_equal(poc_shu2026(model,q)$upper,0)
  extension <- poc_shu2026(model,poc_query(c(x1="y1",x2="y1"),observed_y="y1"))
  expect_match(summary(extension)$bounds$method,"project repeated-outcome extension",fixed=TRUE)
  expect_identical(summary(extension)$bounds$sharpness,"valid_closed_form")
  target <- tempfile("shu-export-")
  poc_export_summary(summary(extension),target)
  exported <- read.delim(file.path(target,"bounds.tsv"),check.names=FALSE)
  expect_identical(exported$method,extension$method)
  expect_identical(exported$sharpness,extension$sharpness)
  expect_error(poc_shu2026(list(),q),"Expected")
  expect_error(poc_shu2026(model,poc_query(c(absent="y1"))),"absent")
  expect_error(poc_shu2026(model,list()),"Expected")
  m <- 30L; n <- 3L
  o <- matrix(1/(m*n),m,n,dimnames=list(paste0("x",1:m),paste0("y",1:n)))
  large <- poc_model(o,o*m,check_lp=FALSE)
  fit <- poc_shu2026(large,poc_query(setNames(rep("y1",m),rownames(o))))
  expect_equal(fit$lower,0)
  expect_true(is.finite(fit$upper))
})
