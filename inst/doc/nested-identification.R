# Executable companion to nested-identification.md. No writes or RNG changes.
library(pnsbib)
observed <- rbind(reference = c(failure = .35, success = .15),
                 active = c(failure = .20, success = .30))
intervention <- rbind(reference = c(failure = .70, success = .30),
                     active = c(failure = .55, success = .45))
query <- poc_query(c(reference = "failure", active = "success"))
model <- poc_model(observed, intervention, lp_backend = "zig")
fits <- list(Li_Pearl = poc_bounds(model, query),
             Shu = poc_shu2026(model, query), LP = poc_exact(model, query))
endpoints <- function(fit) c(lower = fit$lower, upper = fit$upper)
stopifnot(all(vapply(fits, function(fit)
  max(abs(endpoints(fit) - c(.15, .45))) < 1e-8, logical(1))))

# Missing risk under active among factual reference = .2+.1 = .3.
# Missing risk under reference among factual active = .1+.2 = .3.
previous <- c(0, 1)
nested <- lapply(c(1, .5, .25, .1, 0), function(width) {
  constraints <- list(
    poc_constraint(poc_query(c(active = "success"), observed_x = "reference",
                            conditional = TRUE),
                   max(0, .3 - width), min(1, .3 + width)),
    poc_constraint(poc_query(c(reference = "success"), observed_x = "active",
                            conditional = TRUE),
                   max(0, .3 - width), min(1, .3 + width)))
  partial <- poc_partial_model(observed, constraints, lp_backend = "zig")
  result <- endpoints(poc_partial_bounds(partial, query))
  stopifnot(result[1] >= previous[1] - 1e-8,
            result[2] <= previous[2] + 1e-8,
            result[1] <= .3 + 1e-8, .3 <= result[2] + 1e-8)
  previous <<- result
  data.frame(width = width, lower = result[1], upper = result[2])
})
nested <- do.call(rbind, nested)
stopifnot(max(abs(unlist(nested[nrow(nested), c("lower", "upper")]) -
                    c(.15, .45))) < 1e-8)
print(nested, row.names = FALSE)
