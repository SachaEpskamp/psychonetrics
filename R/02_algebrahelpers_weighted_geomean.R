weighted.geomean <- function(x, w, ...) exp(weighted.mean(log(x), w, ...))
geomean <- function(x, ...) exp(mean(log(x), ...))