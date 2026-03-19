# Functions to evaluate the diversity/variation in numeric vectors


# ---- Gini for a numeric vector (nonnegative) ----
gini_vec <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x >= 0]
  n <- length(x); if (n == 0) return(NA_real_)
  sx <- sum(x);  if (sx == 0) return(0)
  x <- sort(x)
  (2 * sum((1:n) * x) / (n * sx)) - (n + 1) / n
}


# ---- Hill index H_m on the top m values of a heavy-tailed vector ----
# Typical choice: m = floor(sqrt(n_nonzero)) or a small fixed m (e.g., 10–50)
hill_index <- function(x, m = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  n <- length(x); if (n < 3) return(NA_real_)
  x <- sort(x, decreasing = TRUE)
  if (is.null(m)) m <- floor(sqrt(n))
  if (m < 1 || m >= n) return(NA_real_)
  k_cut <- x[m + 1]
  1 / mean(log(x[1:m] / k_cut))
}