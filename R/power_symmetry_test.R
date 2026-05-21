#' @import dplyr
#' @import tidyr
#' @import ggpubr
#' @importFrom purrr pmap 
#' @export
#'
#' @title
#' Compute empirical power using Monte Carlo methods
#' @description
#' This function computes the empirical power using Monte Carlo methods for a large number of
#' tests for symmetry with known median. To see the available tests, see details.
#' @return
#' A data frame with empirical powers for each test statistic and distribution.
#' @details
#' Computes empirical power on a tribble or data frame that contains the parameters of the
#' Generalised Lambda Distribution (GLD).
#' Jk: Corzo, J. and Babativa, G. (2013),
#' R: McWiliams, P. (1990),
#' Rs: Baklizi, A. (2003),
#' Mp: Modarres and Gastwirth (1996).
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#' @param data data frame or tribble with the case name and the four lambda parameters of the GLD.
#' @param statis Test statistic to be used. By default \code{statis = c("Bk", "Jk")}.
#'   Available options are \code{"Bk"}, \code{"Bkc"}, \code{"Ck"}, \code{"Jk"}, \code{"R"}, \code{"Rs"}, \code{"Mp"}.
#' @param Bk Integer or vector of integers with the cut-off parameter(s) for the \eqn{B_k} statistic.
#'   By default \code{Bk = 5}. Multiple values can be supplied, e.g. \code{Bk = c(5, 6, 7)}.
#' @param Ck Integer (or vector) with the cut-off(s) for the \eqn{C_k} statistic (conditional Bk).
#'   Requires \eqn{k \geq 10} for reliable size control. By default \code{Ck = 10}.
#' @param Jk Integer or vector of integers with the cut-off parameter(s) for the \eqn{J_k} statistic.
#'   By default \code{Jk = 6}. Multiple values can be supplied, e.g. \code{Jk = c(4, 6, 8)}.
#' @param type Type of test. When the test is with known median, select \code{type = "k"};
#'   use \code{type = "u"} for unknown median (estimated from the data).
#' @param alpha Significance level of the test. By default \eqn{\alpha = 0.05}.
#' @param nsize Sample size for each Monte Carlo replicate.
#' @param rep Number of replicates for the Monte Carlo process.
#' @param plot Logical. If \code{TRUE} (default), draws the empirical power function.
#' @param median Numeric vector with one median value per row in \code{data}. Each element
#'   is the centre parameter used in \code{symmetry_test()} for the corresponding distribution.
#'   By default \code{median = 0} (same value applied to all distributions). You can pass a
#'   column from a pre-computed percentiles table, e.g. \code{median = percentiles$q50}.
#' @param ncores Number of cores to use for parallel computation. By default \code{ncores = NULL},
#'   which automatically uses all available cores minus one.
#' @references
#' Corzo, J., & Babativa, G. (2013). A modified runs test for symmetry. Journal of Statistical
#' Computation and Simulation, 83(5), 984-991.
#' @examples
#' \dontrun{
#' distributions <- tibble::tribble(
#'   ~case,     ~lambda1, ~lambda2, ~lambda3, ~lambda4,
#'   "Case 1A",        0, 0.197454, 0.134915, 0.134915
#' )
#' P20 <- power_symmetry_test(data = distributions, statis = c("Bk", "Jk"),
#'                            Bk = 5, Jk = 6, alpha = 0.05, nsize = 20, rep = 1000,
#'                            plot = FALSE, median = 0)
#'
#' # Using per-distribution medians from a percentiles table:
#' P30 <- power_symmetry_test(data   = distributions,
#'                            statis = c("Bk", "Jk", "R", "Rs", "Mp"),
#'                            Bk     = c(5, 10),
#'                            alpha  = 0.05,
#'                            nsize  = 30,
#'                            rep    = 3000,
#'                            plot   = FALSE,
#'                            median = percentiles$q50,
#'                            ncores = parallel::detectCores() - 1)
#' }

power_symmetry_test <- function(data, statis = c("Bk", "Jk", "R", "Rs", "Mp"),
                                Bk = 5, Ck = 10, Jk = 6, type = "k",
                                alpha = 0.05, nsize, rep, plot = TRUE,
                                median = 0,
                                ncores = NULL) {

  ni <- nsize
  B  <- rep
  n_dist <- nrow(data)

  # --- Expandir median a un vector de longitud n_dist ---
  if (length(median) == 1L) {
    median_vec <- rep(median, n_dist)
  } else if (length(median) == n_dist) {
    median_vec <- as.numeric(median)
  } else {
    stop("`median` must be either a single value or a numeric vector with one ",
         "element per row in `data` (", n_dist, " rows).")
  }

  # --- Número de cores automático si no se especifica ---
  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)

  # --- Precalcular valores críticos UNA sola vez (fuera de las réplicas) ---
  crit <- list()
  if ("Bk" %in% statis) {
    for (k in Bk) {
      qb  <- stats::qbinom(alpha, k, 0.5)
      sum_antes <- if (qb == 0) 0 else sum(stats::dbinom(0:(qb - 1), k, 0.5))
      acc <- (alpha - sum_antes) / stats::dbinom(qb, k, 0.5)
      crit[[paste0("Bk", k)]] <- list(q = qb, acc = acc)
    }
  }
  if ("R" %in% statis) {
    qr  <- stats::qbinom(alpha, ni - 1, 0.5)
    acc <- (alpha - sum(stats::dbinom(0:(qr - 1), ni - 1, 0.5))) /
      stats::dbinom(qr, ni - 1, 0.5)
    crit[["R"]] <- list(q = qr, acc = acc)
  }

  # --- Función de rechazo vectorizada ---

  compute_reject <- function(stat_name, stat_value, p_value, n1 = NA, n0 = NA) {
    sv <- stat_value - 1

    # Bk — valor crítico binomial precalculado, test aleatorizado
    if (grepl("^Bk[0-9]", stat_name)) {
      info <- crit[[stat_name]]
      return(ifelse(sv < info$q, 1,
                    ifelse(sv == info$q, info$acc, 0)))
    }
    # R — valor crítico binomial precalculado, test aleatorizado
    if (stat_name == "R") {
      info <- crit[["R"]]
      return(ifelse(sv < info$q, 1,
                    ifelse(sv == info$q, info$acc, 0)))
    }
    # Bkc y Ck — aleatorización condicional muestra a muestra
    if (stat_name == "Bkc" || grepl("^Ck[0-9]", stat_name)) {
      return(mapply(function(rv, n1i, n0i) {
        if (is.na(n1i) || is.na(n0i) || n1i == 0L || n0i == 0L) return(0)
        r_vals <- seq_len(n1i + n0i)
        probs  <- randtests::druns(r_vals - 1, n1i, n0i)
        cdf    <- cumsum(probs)
        qi_idx <- suppressWarnings(max(which(cdf <= alpha)))
        if (!is.finite(qi_idx) || qi_idx == 0) return(0)
        qi   <- r_vals[qi_idx]
        acci <- (alpha - cdf[qi_idx]) / probs[qi_idx + 1]
        ifelse(rv < qi, 1, ifelse(rv == qi, acci, 0))
      }, stat_value, n1, n0))
    }
    # Jk*, Rs, Mp* — rechazo directo por p-valor
    return(as.integer(p_value < alpha))
  }

  lsamples <- purrr::pmap(data[, 2:5], samples, ni = ni, B = B)

  # --- Clúster paralelo ---
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl,
                          varlist = c("symmetry_test", "statis", "Bk", "Jk", "type",
                                      "alpha", "ni", "B", "crit", "compute_reject",
                                      "median_vec"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(randtests)
  })

  # --- Paralelizado por distribución ---
  # Se pasa el índice i junto con las muestras para recuperar la mediana correcta
  outPower <- parallel::parLapply(cl, seq_along(lsamples), function(i) {

    x   <- lsamples[[i]]
    med <- median_vec[[i]]   # mediana específica de esta distribución

    result <- apply(x, 2, function(y) {
      symmetry_test(y, statis = statis, Bk = Bk, Ck = Ck, Jk = Jk, type = type, median = med)
    })

    dplyr::bind_rows(result) %>%
      group_by(stat) %>%
      mutate(reject = compute_reject(stat[1], stat.value, p.value, n1, n0)) %>%
      summarise(power = mean(reject), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = stat, values_from = power)
  })

  Power <- dplyr::bind_rows(outPower) %>%
    cbind(data[, "case"])

  if (plot) print(Power.plot(Power, ni))

  return(Power)
}
