#' Simulate a Random-Intercept-and-Slope GLMM
#'
#' Simulates clustered responses from a negative-binomial or compound
#' Poisson–Gamma GLMM with log link and Gaussian random effects.
#' Singular positive-semidefinite random-effects covariance matrices are
#' supported. Uses the current R random-number generator state; call
#' `set.seed()` externally for reproducibility.
#'
#' @param n_clusters Positive integer number of clusters.
#' @param cluster_size Integer observations per cluster (at least two).
#' @param beta Numeric vector of length two (intercept and slope).
#' @param Gamma Symmetric positive-semidefinite 2-by-2 covariance matrix.
#' @param family Either `"nbinom"` or `"cpg"`.
#' @param family_params Named numeric vector: `theta` for negative binomial,
#'   or `phi` and `power` for compound Poisson–Gamma.
#' @param x Optional finite numeric vector of length `n_clusters * cluster_size`.
#'   If supplied, overrides `x_design`.
#' @param x_design Either `"binary_balanced"` or `"continuous"`.
#'
#' @return A data frame with factor `id`, numeric `x`, and response `y`.
#' @export
#' @examples
#' set.seed(123)
#' d <- simulate_glmm(4, 6, c(0.2, 0.3), diag(c(0.1, 0.2)),
#'                    family = "nbinom", family_params = c(theta = 2))
#' head(d)

simulate_glmm <- function(
    n_clusters,
    cluster_size,
    beta,
    Gamma,
    family = c("nbinom", "cpg"),
    family_params,
    x = NULL,
    x_design = c("binary_balanced", "continuous")
) {
  
  family <- match.arg(family)
  x_design <- match.arg(x_design)
  
  # ----------------------------------------------------------
  # 1. Validate dimensions and parameters
  # ----------------------------------------------------------
  
  stopifnot(
    length(n_clusters) == 1L,
    is.numeric(n_clusters),
    is.finite(n_clusters),
    n_clusters >= 1,
    n_clusters == floor(n_clusters),
    
    length(cluster_size) == 1L,
    is.numeric(cluster_size),
    is.finite(cluster_size),
    cluster_size >= 2,
    cluster_size == floor(cluster_size),
    
    is.numeric(beta),
    length(beta) == 2L,
    all(is.finite(beta)),
    
    is.matrix(Gamma),
    is.numeric(Gamma),
    identical(dim(Gamma), c(2L, 2L)),
    all(is.finite(Gamma))
  )
  
  if (!isTRUE(all.equal(
    Gamma, t(Gamma), tolerance = 1e-10
  ))) {
    stop("Gamma must be symmetric.")
  }
  
  Gamma <- (Gamma + t(Gamma)) / 2
  
  eig <- eigen(Gamma, symmetric = TRUE)
  
  if (min(eig$values) < -1e-10) {
    stop("Gamma must be positive semidefinite.")
  }
  
  # ----------------------------------------------------------
  # 2. Validate family-specific parameters
  # ----------------------------------------------------------
  
  pn <- switch(
    family,
    nbinom = "theta",
    cpg = c("phi", "power")
  )
  
  if (!is.numeric(family_params) ||
      is.null(names(family_params)) ||
      length(family_params) != length(pn) ||
      anyDuplicated(names(family_params)) ||
      !setequal(names(family_params), pn) ||
      any(!is.finite(family_params))) {
    
    stop(
      "family_params must contain exactly: ",
      paste(pn, collapse = ", ")
    )
  }
  
  family_params <- family_params[pn]
  
  if (family == "nbinom") {
    
    if (family_params[["theta"]] <= 0) {
      stop("Negative binomial requires theta > 0.")
    }
    
  } else if (family == "cpg") {
    
    if (family_params[["phi"]] <= 0 ||
        family_params[["power"]] <= 1 ||
        family_params[["power"]] >= 2) {
      stop("CPG requires phi > 0 and 1 < power < 2.")
    }
  }
  
  # ----------------------------------------------------------
  # 3. Cluster identifiers
  # ----------------------------------------------------------
  
  n <- n_clusters * cluster_size
  
  if (!is.finite(n) || n > .Machine$integer.max) {
    stop("Total number of observations is too large.")
  }
  
  n <- as.integer(n)
  n_clusters <- as.integer(n_clusters)
  cluster_size <- as.integer(cluster_size)
  
  id <- rep(
    seq_len(n_clusters),
    each = cluster_size
  )
  
  # ----------------------------------------------------------
  # 4. Generate covariate design
  # ----------------------------------------------------------
  
  if (is.null(x)) {
    
    if (x_design == "binary_balanced") {
      
      n0 <- floor(cluster_size / 2)
      n1 <- cluster_size - n0
      
      x <- unlist(
        lapply(seq_len(n_clusters), function(i) {
          
          if (cluster_size %% 2L == 1L &&
              sample.int(2L, 1L) == 1L) {
            counts <- c(n1, n0)
          } else {
            counts <- c(n0, n1)
          }
          
          sample(
            rep(c(0, 1), times = counts),
            size = cluster_size,
            replace = FALSE
          )
        }),
        use.names = FALSE
      )
      
    } else if (x_design == "continuous") {
      
      x <- rep(
        seq(0, 1, length.out = cluster_size),
        times = n_clusters
      )
      
    }
  }
  
  if (!is.numeric(x) ||
      length(x) != n ||
      any(!is.finite(x))) {
    
    stop(
      "x must be a finite numeric vector of length n."
    )
  }
  
  # ----------------------------------------------------------
  # 5. Generate random intercepts and slopes
  # ----------------------------------------------------------
  
  A <- eig$vectors %*%
    diag(
      sqrt(pmax(eig$values, 0)),
      nrow = 2L
    )
  
  B <- matrix(
    stats::rnorm(n_clusters * 2L),
    nrow = n_clusters,
    ncol = 2L
  ) %*% t(A)
  
  # ----------------------------------------------------------
  # 6. Linear predictor and conditional mean
  # ----------------------------------------------------------
  
  eta <- beta[1L] +
    beta[2L] * x +
    B[id, 1L] +
    B[id, 2L] * x
  
  mu <- exp(eta)
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("Invalid conditional means.")
  }
  
  # ----------------------------------------------------------
  # 7. Generate family-specific responses
  # ----------------------------------------------------------
  
  y <- switch(
    family,
    
    nbinom = {
      
      stats::rnbinom(
        n = n,
        size = family_params[["theta"]],
        mu = mu
      )
      
    },
    
    cpg = {
      
      phi <- family_params[["phi"]]
      power <- family_params[["power"]]
      
      lambda <- mu^(2 - power) /
        (phi * (2 - power))
      
      alpha <- (2 - power) /
        (power - 1)
      
      scale <- phi * (power - 1) *
        mu^(power - 1)
      
      N <- stats::rpois(n, lambda = lambda)
      
      y <- numeric(n)
      
      positive <- N > 0L
      
      y[positive] <- stats::rgamma(
        n = sum(positive),
        shape = N[positive] * alpha,
        scale = scale[positive]
      )
      
      y
    }
  )
  
  # ----------------------------------------------------------
  # 8. Return simulated data
  # ----------------------------------------------------------
  
  base::data.frame(
    id = factor(id),
    x = as.numeric(x),
    y = y
  )
}
