
# ============================================================
# DETERMINISTIC FISHER INFORMATION FOR GLMMs
# At the covariance boundary Gamma = 0
#
# Supported families:
#   1. Negative binomial, log link
#   2. Compound Poisson-Gamma, log link
#
# Public function:
#   fim_glmm()
#
# Family-specific parameters:
#   NB:  c(theta = ...)
#   CPG: c(phi = ..., power = ...)
#
# Dependencies for CPG:
#   tweedie
#
# No Monte Carlo simulation is used.
# ============================================================


# ============================================================
# 1. PARAMETER VALIDATION
# ============================================================

nb_validate_params <- function(theta) {
  
  if (!is.numeric(theta) ||
      length(theta) != 1L ||
      !is.finite(theta) ||
      theta <= 0) {
    stop("'theta' must be a finite positive scalar.")
  }
  
  invisible(TRUE)
}


cp_validate_params <- function(phi, power) {
  
  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("'phi' must be a finite positive scalar.")
  }
  
  if (!is.numeric(power) ||
      length(power) != 1L ||
      !is.finite(power) ||
      power <= 1 ||
      power >= 2) {
    stop("'power' must satisfy 1 < power < 2.")
  }
  
  invisible(TRUE)
}


# ============================================================
# 2. NEGATIVE BINOMIAL NUISANCE INFORMATION
# ============================================================

.nb_theta_info <- function(
    mu,
    theta,
    tail_tol = 1e-12,
    max_terms = 100000L
) {
  
  nb_validate_params(theta)
  
  K <- min(64L, max_terms)
  
  repeat {
    
    k <- seq.int(0L, K - 1L)
    
    surv <- stats::pnbinom(
      k,
      size = theta,
      mu = mu,
      lower.tail = FALSE
    )
    
    value <- sum(
      surv / (theta + k)^2
    ) - 1 / theta + 1 / (theta + mu)
    
    bound <- surv[K] * trigamma(theta + K)
    
    if (bound <= tail_tol) {
      return(value)
    }
    
    if (K >= max_terms) {
      stop("NB theta-information series failed to converge.")
    }
    
    K <- min(max_terms, 2L * K)
  }
}


# ============================================================
# 3. CPG OBSERVATION-LEVEL NUISANCE SCORES
# ============================================================

# For y > 0, the CPG density is an infinite
# Poisson-Gamma mixture.
#
# We evaluate its phi and power scores using
# normalized mixture weights.

cp_nuisance_score <- function(
    y,
    mu,
    phi,
    power,
    tol = 1e-12,
    max_terms = 100000L
) {
  
  cp_validate_params(phi, power)
  
  if (!is.numeric(y) ||
      !is.numeric(mu) ||
      any(!is.finite(y)) ||
      any(!is.finite(mu)) ||
      any(y < 0) ||
      any(mu <= 0)) {
    stop("Invalid y or mu.")
  }
  
  if (length(mu) == 1L) {
    mu <- rep(mu, length(y))
  }
  
  if (length(y) != length(mu)) {
    stop("y and mu must have compatible lengths.")
  }
  
  out <- matrix(
    0,
    nrow = length(y),
    ncol = 2L,
    dimnames = list(NULL, c("phi", "power"))
  )
  
  for (i in seq_along(y)) {
    
    yi <- y[i]
    mui <- mu[i]
    
    lambda <- mui^(2 - power) /
      (phi * (2 - power))
    
    alpha <- (2 - power) / (power - 1)
    
    scale <- phi * (power - 1) *
      mui^(power - 1)
    
    Lp <- -log(mui) + 1 / (2 - power)
    
    Sp <- log(mui) + 1 / (power - 1)
    
    alpha_prime <- -1 / (power - 1)^2
    
    # Exact atom at zero
    if (yi == 0) {
      
      out[i, ] <- c(
        phi = lambda / phi,
        power = -lambda * Lp
      )
      
      next
    }
    
    # Positive observations
    log_y <- log(yi)
    log_lambda <- log(lambda)
    log_scale <- log(scale)
    
    logw <- numeric(0)
    score_phi <- numeric(0)
    score_power <- numeric(0)
    
    max_logw <- -Inf
    low_count <- 0L
    
    for (n in seq_len(max_terms)) {
      
      na <- n * alpha
      
      lw <- n * log_lambda -
        lgamma(n + 1) -
        lgamma(na) +
        (na - 1) * log_y -
        na * log_scale
      
      uphi <- (
        lambda - n - na + yi / scale
      ) / phi
      
      upower <- (
        (n - lambda) * Lp +
          n * alpha_prime *
          (log_y - log_scale - digamma(na)) +
          (yi / scale - na) * Sp
      )
      
      logw[n] <- lw
      score_phi[n] <- uphi
      score_power[n] <- upower
      
      max_logw <- max(max_logw, lw)
      
      if (lw < max_logw + log(tol)) {
        low_count <- low_count + 1L
      } else {
        low_count <- 0L
      }
      
      if (low_count >= 30L && n > 30L) {
        break
      }
      
      if (n == max_terms) {
        stop("CPG mixture series failed to converge.")
      }
    }
    
    weights <- exp(logw - max_logw)
    weights <- weights / sum(weights)
    
    out[i, 1] <- sum(weights * score_phi)
    out[i, 2] <- sum(weights * score_power)
  }
  
  out
}


# ============================================================
# 4. CPG NUISANCE FISHER INFORMATION
# ============================================================

# Computes E[u_psi u_psi^T] deterministically.
#
# The atom at zero is handled exactly.
# The positive component uses numerical quadrature.
#
# Note:
# Finite integration limits and mixture truncation
# require convergence checks for research use.

.cp_psi_info <- function(
    mu,
    phi,
    power,
    rel_tol = 1e-7,
    tail_prob = 1e-12,
    max_terms = 100000L,
    subdivisions = 300L
) {
  
  cp_validate_params(phi, power)
  
  if (!requireNamespace("tweedie", quietly = TRUE)) {
    stop("Package 'tweedie' is required for CPG information.")
  }
  
  lambda <- mu^(2 - power) /
    (phi * (2 - power))
  
  p0 <- exp(-lambda)
  
  u0 <- c(
    phi = lambda / phi,
    power = -lambda *
      (-log(mu) + 1 / (2 - power))
  )
  
  I <- p0 * tcrossprod(u0)
  
  # Choose integration limits using the positive
  # component of the Tweedie distribution.
  # tail_prob <- min(1e-10, rel_tol * 1e-3)
  
  lo <- tweedie::qtweedie(
    p = p0 + tail_prob * (1 - p0),
    mu = mu,
    phi = phi,
    power = power
  )
  
  hi <- tweedie::qtweedie(
    p = 1 - tail_prob,
    mu = mu,
    phi = phi,
    power = power
  )
  
  if (!is.finite(lo) ||
      !is.finite(hi) ||
      lo <= 0 ||
      hi <= lo) {
    stop("Unable to establish CPG integration bounds.")
  }
  
  integrand <- function(t, component) {
    
    vapply(t, function(tt) {
      
      y <- exp(tt)
      
      if (!is.finite(y) || y <= 0) {
        return(0)
      }
      
      density <- tweedie::dtweedie(
        y = y,
        mu = mu,
        phi = phi,
        power = power
      )
      
      if (!is.finite(density)) {
        stop("Non-finite CPG density.")
      }
      
      if (density == 0) {
        return(0)
      }
      
      u <- cp_nuisance_score(
        y = y,
        mu = mu,
        phi = phi,
        power = power,
        max_terms = max_terms
      )[1, ]
      
      value <- switch(
        component,
        "11" = u[1]^2,
        "12" = u[1] * u[2],
        "22" = u[2]^2
      )
      
      value * density * y
    }, numeric(1))
  }
  
  components <- c("11", "12", "22")
  
  values <- vapply(
    components,
    function(component) {
      
      stats::integrate(
        f = function(t) integrand(t, component),
        lower = log(lo),
        upper = log(hi),
        rel.tol = rel_tol,
        subdivisions = subdivisions
      )$value
    },
    numeric(1)
  )
  
  I[1, 1] <- I[1, 1] + values["11"]
  I[1, 2] <- I[1, 2] + values["12"]
  I[2, 1] <- I[1, 2]
  I[2, 2] <- I[2, 2] + values["22"]
  
  dimnames(I) <- list(
    c("phi", "power"),
    c("phi", "power")
  )
  
  I
}


# ============================================================
# 5. FAMILY SPECIFICATIONS
# ============================================================


nbinom_theoretical <- function() {
  
  list(
    name = "nbinom",
    param_names = "theta",
    
    validate = function(params) {
      nb_validate_params(params[["theta"]])
    },
    
    moments = function(mu, params) {
      
      theta <- as.numeric(params[["theta"]])
      mu <- as.numeric(mu)
      
      v <- theta * mu / (theta + mu)
      
      C <- theta * mu * (theta + 3 * mu) /
        (theta + mu)^2
      
      D <- stats::setNames(
        as.numeric(-mu^2 / (theta + mu)^2),
        "theta"
      )
      
      list(
        A = as.numeric(v),
        B = as.numeric(v),
        C = as.numeric(C),
        D = D
      )
    },
    
    nuisance_info = function(
    mu, params, tail_tol, max_terms
    ) {
      
      matrix(
        .nb_theta_info(
          mu = mu,
          theta = unname(params[["theta"]]),
          tail_tol = tail_tol,
          max_terms = max_terms
        ),
        nrow = 1L,
        dimnames = list("theta", "theta")
      )
    }
  )
}

cpg_theoretical <- function() {
  
  list(
    name = "cpg",
    param_names = c("phi", "power"),
    
    validate = function(params) {
      
      cp_validate_params(
        phi = params[["phi"]],
        power = params[["power"]]
      )
    },
    
    moments = function(mu, params) {
      
      phi <- as.numeric(params[["phi"]])
      power <- as.numeric(params[["power"]])
      mu <- as.numeric(mu)
      
      v <- mu^(2 - power) / phi
      
      D <- stats::setNames(
        as.numeric(c(
          v / phi,
          v * log(mu)
        )),
        c("phi", "power")
      )
      
      list(
        A = as.numeric(v),
        B = as.numeric(v),
        C = as.numeric((power^2 - power + 1) * v),
        D = D
      )
    },
    
    nuisance_info = function(
    mu, params, tail_tol, max_terms
    ) {
      
      .cp_psi_info(
        mu = mu,
        phi = params[["phi"]],
        power = params[["power"]],
        rel_tol = tail_tol,
        max_terms = max_terms
      )
    }
  )
}


theoretical_family <- function(family) {
  
  if (is.list(family) &&
      is.function(family$moments) &&
      is.function(family$nuisance_info) &&
      is.function(family$validate)) {
    return(family)
  }
  
  if (!is.character(family) ||
      length(family) != 1L ||
      is.na(family)) {
    stop("Invalid theoretical family.")
  }
  
  switch(
    family,
    nbinom = nbinom_theoretical(),
    nbinom_log = nbinom_theoretical(),
    cpg = cpg_theoretical(),
    cp = cpg_theoretical(),
    stop("Unsupported theoretical family: ", family)
  )
}


# ============================================================
# 6. INTERNAL FISHER INFORMATION ENGINE
# ============================================================

.fim_glmm_matrix <- function(
    X,
    Z,
    beta,
    family,
    family_params,
    tail_tol = 1e-8,
    max_terms = 100000L
) {
  
  spec <- theoretical_family(family)
  
  pn <- spec$param_names
  k <- length(pn)
  
  if (!is.numeric(family_params) ||
      is.null(names(family_params)) ||
      length(family_params) != k ||
      anyDuplicated(names(family_params)) ||
      !setequal(names(family_params), pn) ||
      any(!is.finite(family_params))) {
    stop("Invalid family_params.")
  }
  
  family_params <- family_params[pn]
  spec$validate(family_params)
  
  if (!is.list(X) || !is.list(Z) ||
      length(X) == 0L ||
      length(X) != length(Z)) {
    stop("X and Z must be nonempty equal-length lists.")
  }
  
  if (!is.matrix(X[[1L]]) ||
      !is.matrix(Z[[1L]])) {
    stop("X and Z must contain matrices.")
  }
  
  p <- ncol(X[[1L]])
  q <- ncol(Z[[1L]])
  
  if (p < 1L || q < 1L ||
      !is.numeric(beta) ||
      length(beta) != p ||
      any(!is.finite(beta))) {
    stop("Invalid beta or design dimensions.")
  }
  
  pairs <- which(
    upper.tri(matrix(0, q, q), diag = TRUE),
    arr.ind = TRUE
  )
  
  ng <- nrow(pairs)
  
  factors <- ifelse(
    pairs[, 1L] == pairs[, 2L],
    0.5, 1
  )
  
  bn <- colnames(X[[1L]])
  
  if (is.null(bn)) {
    bn <- paste0("b", seq_len(p) - 1L)
  }
  
  gn <- paste0(
    "g", pairs[, 1L], pairs[, 2L]
  )
  
  nm <- c(bn, gn, pn)
  
  bi <- seq_len(p)
  gi <- p + seq_len(ng)
  fi <- p + ng + seq_len(k)
  
  I <- matrix(
    0,
    nrow = length(nm),
    ncol = length(nm),
    dimnames = list(nm, nm)
  )
  
  cache <- new.env(parent = emptyenv())
  
  for (i in seq_along(X)) {
    
    Xi <- X[[i]]
    Zi <- Z[[i]]
    
    if (!is.matrix(Xi) ||
        !is.matrix(Zi) ||
        ncol(Xi) != p ||
        ncol(Zi) != q ||
        nrow(Xi) != nrow(Zi) ||
        nrow(Xi) < 1L ||
        any(!is.finite(Xi)) ||
        any(!is.finite(Zi))) {
      stop("Invalid matrices in cluster ", i)
    }
    
    mu <- drop(exp(Xi %*% beta))
    
    if (any(!is.finite(mu)) ||
        any(mu <= 0)) {
      stop("Invalid mean in cluster ", i)
    }
    
    n <- length(mu)
    
    A <- numeric(n)
    B <- numeric(n)
    C <- numeric(n)
    
    D <- matrix(
      0,
      nrow = n,
      ncol = k,
      dimnames = list(NULL, pn)
    )
    
    for (j in seq_len(n)) {
      
      moments <- spec$moments(
        mu[j], family_params
      )
      
      A[j] <- moments$A
      B[j] <- moments$B
      C[j] <- moments$C
      
      if (length(moments$D) != k ||
          is.null(names(moments$D)) ||
          !setequal(names(moments$D), pn)) {
        stop("Invalid family cross-information.")
      }
      
      D[j, ] <- moments$D[pn]
      
      key <- sprintf("%.17g", mu[j])
      
      if (!exists(
        key, envir = cache,
        inherits = FALSE
      )) {
        
        F <- spec$nuisance_info(
          mu = mu[j],
          params = family_params,
          tail_tol = tail_tol,
          max_terms = max_terms
        )
        
        if (!is.matrix(F) ||
            !identical(dim(F), c(k, k)) ||
            any(!is.finite(F))) {
          stop("Invalid family nuisance information.")
        }
        
        assign(key, F, envir = cache)
      }
      
      I[fi, fi] <- I[fi, fi] +
        get(key, envir = cache)
    }
    
    # Beta-beta information
    I[bi, bi] <- I[bi, bi] +
      crossprod(Xi, A * Xi)
    
    # Random-effect score covariance
    M <- crossprod(Zi, A * Zi)
    
    for (u in seq_len(ng)) {
      
      a <- pairs[u, 1L]
      b <- pairs[u, 2L]
      
      zab <- Zi[, a] * Zi[, b]
      
      # Beta-gamma information
      I[bi, gi[u]] <- I[bi, gi[u]] +
        factors[u] *
        drop(crossprod(Xi, B * zab))
      
      # Gamma-family information
      I[gi[u], fi] <- I[gi[u], fi] +
        factors[u] *
        colSums(D * zab)
      
      # Gamma-gamma information
      for (h in seq_len(ng)) {
        
        c <- pairs[h, 1L]
        d <- pairs[h, 2L]
        
        I[gi[u], gi[h]] <-
          I[gi[u], gi[h]] +
          factors[u] * factors[h] * (
            M[a, c] * M[b, d] +
              M[a, d] * M[b, c] +
              sum(
                C * zab *
                  Zi[, c] * Zi[, d]
              )
          )
      }
    }
  }
  
  # Complete symmetric blocks
  I[gi, bi] <- t(
    I[bi, gi, drop = FALSE]
  )
  
  I[fi, gi] <- t(
    I[gi, fi, drop = FALSE]
  )
  
  # Beta-family information is zero
  # for the supported mean parameterizations.
  
  # Efficient information for Gamma
  nuisance <- c(bi, fi)
  
  Inuis <- I[
    nuisance, nuisance,
    drop = FALSE
  ]
  
  Ig_nuis <- I[
    gi, nuisance,
    drop = FALSE
  ]
  
  metric <- I[gi, gi, drop = FALSE] -
    Ig_nuis %*%
    solve(Inuis, t(Ig_nuis))
  
  metric <- (metric + t(metric)) / 2
  
  list(
    I_full = I,
    I_metric = metric,
    beta_idx = bi,
    gamma_idx = gi,
    family_idx = fi,
    family_param_names = pn,
    method = "theoretical",
    boundary = TRUE
  )
}


# ============================================================
# 7. PUBLIC FORMULA-BASED INTERFACE
# ============================================================

#' Fisher information for a GLMM at zero random-effects covariance
#'
#' Computes deterministic expected Fisher information at the boundary
#' `Gamma = 0` for supported log-link GLMM families. Returns the full
#' information and the efficient covariance-parameter information obtained
#' by eliminating fixed-effect and family-specific nuisance parameters.
#'
#' @param fixed One-sided fixed-effects model formula.
#' @param random One-sided random-effects model formula.
#' @param data Nonempty design data frame.
#' @param cluster Name of the cluster identifier column in `data`.
#' @param beta Numeric fixed-effect coefficient vector in model-matrix order.
#' @param family Either `"nbinom"` or `"cpg"` (supported aliases accepted),
#'   or an internal family specification.
#' @param family_params Named numeric vector: `c(theta = ...)` for negative
#'   binomial or `c(phi = ..., power = ...)` for compound Poisson-Gamma.
#' @param tail_tol Positive numerical tolerance for nuisance information.
#' @param max_terms Positive maximum series length for numerical summation.
#'
#' @return A list with `I_full`, `I_metric`, parameter block indices,
#'   family metadata, and design metadata. `I_metric` is the Schur
#'   complement of the nuisance-information block.
#' @details The CPG nuisance-information integration and truncated mixture
#'   require numerical convergence checks. This function does not estimate
#'   parameters or use Monte Carlo simulation.
#' @examples
#' design <- expand.grid(x = c(0, 1), id = seq_len(5))
#' info <- fim_glmm(~ x, ~ 1 + x, design, "id", c(0.5, 0.2),
#'                  "nbinom", c(theta = 2))
#' @export
fim_glmm <- function(
    fixed,
    random,
    data,
    cluster,
    beta,
    family,
    family_params,
    tail_tol = 1e-8,
    max_terms = 100000L
) {
  
  # Validate formulas
  if (!inherits(fixed, "formula") ||
      length(fixed) != 2L) {
    stop("'fixed' must be a one-sided formula.")
  }
  
  if (!inherits(random, "formula") ||
      length(random) != 2L) {
    stop("'random' must be a one-sided formula.")
  }
  
  # Validate data
  if (!is.data.frame(data) ||
      nrow(data) == 0L) {
    stop("'data' must be a nonempty data frame.")
  }
  
  if (!is.character(cluster) ||
      length(cluster) != 1L ||
      is.na(cluster) ||
      !cluster %in% names(data)) {
    stop("'cluster' must name a column in 'data'.")
  }
  
  if (anyNA(data[[cluster]])) {
    stop("Cluster identifiers cannot contain missing values.")
  }
  
  if (!is.numeric(tail_tol) ||
      length(tail_tol) != 1L ||
      !is.finite(tail_tol) ||
      tail_tol <= 0 ||
      tail_tol >= 1) {
    stop("'tail_tol' must be between 0 and 1.")
  }
  
  if (!is.numeric(max_terms) ||
      length(max_terms) != 1L ||
      !is.finite(max_terms) ||
      max_terms < 1 ||
      max_terms != as.integer(max_terms)) {
    stop("'max_terms' must be a positive integer.")
  }
  
  max_terms <- as.integer(max_terms)
  
  # Construct design matrices
  X <- stats::model.matrix(
    object = fixed,
    data = data,
    na.action = na.fail
  )
  
  Z <- stats::model.matrix(
    object = random,
    data = data,
    na.action = na.fail
  )
  
  # Validate beta
  if (!is.numeric(beta) ||
      length(beta) != ncol(X) ||
      any(!is.finite(beta))) {
    stop(
      "'beta' must be a finite numeric vector of length ",
      ncol(X), "."
    )
  }
  
  if (ncol(Z) == 0L) {
    stop("The random-effects design must be nonempty.")
  }
  
  # Resolve family
  spec <- theoretical_family(family)
  
  pn <- spec$param_names
  
  if (!is.numeric(family_params) ||
      is.null(names(family_params)) ||
      length(family_params) != length(pn) ||
      anyDuplicated(names(family_params)) ||
      !setequal(names(family_params), pn) ||
      any(!is.finite(family_params))) {
    stop(
      "'family_params' must contain exactly: ",
      paste(pn, collapse = ", ")
    )
  }
  
  family_params <- family_params[pn]
  
  # Unified family validation
  spec$validate(family_params)
  
  # Split design matrices by cluster
  cluster_id <- data[[cluster]]
  
  idx <- split(
    seq_len(nrow(data)),
    factor(
      cluster_id,
      levels = unique(cluster_id)
    ),
    drop = TRUE
  )
  
  X_list <- lapply(
    idx,
    function(i) X[i, , drop = FALSE]
  )
  
  Z_list <- lapply(
    idx,
    function(i) Z[i, , drop = FALSE]
  )
  
  # Compute Fisher information
  result <- .fim_glmm_matrix(
    X = X_list,
    Z = Z_list,
    beta = beta,
    family = spec,
    family_params = family_params,
    tail_tol = tail_tol,
    max_terms = max_terms
  )
  
  # Metadata
  result$family <- spec$name
  result$family_params <- family_params
  result$fixed <- fixed
  result$random <- random
  result$cluster <- cluster
  result$n_clusters <- length(idx)
  result$n_observations <- nrow(data)
  
  result
}


# ============================================================
# 8. EXAMPLE USAGE
# ============================================================

# design <- expand.grid(
#   x = seq(0, 1, length.out = 5),
#   id = seq_len(20)
# )
#
# fim_nb <- fim_glmm(
#   fixed = ~ x,
#   random = ~ 1 + x,
#   data = design,
#   cluster = "id",
#   beta = c(0.5, 0.2),
#   family = "nbinom",
#   family_params = c(theta = 2)
# )
#
# fim_cpg <- fim_glmm(
#   fixed = ~ x,
#   random = ~ 1 + x,
#   data = design,
#   cluster = "id",
#   beta = c(0.5, 0.2),
#   family = "cpg",
#   family_params = c(phi = 1.2, power = 1.5),
#   tail_tol = 1e-7
# )
#
# fim_nb$I_full
# fim_nb$I_metric
#
# fim_cpg$I_full
# fim_cpg$I_metric
