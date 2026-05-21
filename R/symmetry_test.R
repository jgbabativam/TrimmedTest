#' @import dplyr
#' @importFrom randtests druns
#' @export
#'
#' @title
#' Perform symmetry test
#' @description
#' This function contains a large number of tests for symmetry with known median. To see the available tests, see
#' details.
#' @return
#' A data frame with statistical values and p.values
#' @details
#' Perform symmetry tests on numeric vectors or objects, choose test with vector stat.
#' Jk: Corzo, J. and Babativa, G. (2013),
#' R: McWiliams, P. (1990),
#' Rs: Baklizi, A. (2003)
#' Mp: Modarres and Gastwirth (1996).
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#' @param x vector of numeric information.
#' @param statis Test statistic to be used. By default \code{statis = c("Bk", "Jk")}.
#'   Available options are \code{"Bk"}, \code{"Ck"}, \code{"Bkc"}, \code{"Jk"}, \code{"R"}, \code{"Rs"}, \code{"Mp"}.
#' @param Bk Integer or vector of integers with the cut-off parameter(s) for the \eqn{B_k} statistic.
#'   By default \code{Bk = 5}. Multiple values can be supplied, e.g. \code{Bk = c(5, 6, 7)}.
#' @param Jk Integer or vector of integers with the cut-off parameter(s) for the \eqn{J_k} statistic.
#'   By default \code{Jk = 6}. Multiple values can be supplied, e.g. \code{Jk = c(4, 6, 8)}.
#' @param Ck Integer (or vector) with the cut-off(s) for the \eqn{C_k} statistic.
#'   \eqn{C_k} is the conditional version of \eqn{B_k}: it uses the runs distribution
#'   conditioned on \eqn{(n_{1k}, n_{0k})} among the \eqn{k} most extreme observations,
#'   making it robust to misspecified medians. Requires \eqn{k \\geq 10} for reliable
#'   size control. By default \code{Ck = 10}.
#' @param Bkc Integer with the cut-off parameter for the \eqn{B_{kc}} statistic, which uses
#'   the conditional runs distribution. By default \code{Bkc = 11}.
#' @param Mp Numeric vector with percentile trimming values (in percent) for the \eqn{M_p} statistic.
#'   By default \code{Mp = c(10, 20, 25)}.
#' @param median The centre parameter around which to test symmetry. By default \code{median = 0}.
#' @param type Type of test. When the test is with known median, select \code{type = "k"};
#'   use \code{type = "u"} for unknown median (estimated from the data).
#' @references
#' Corzo, J., & Babativa, G. (2013). A modified runs test for symmetry. Journal of Statistical Computation and Simulation, 83(5), 984-991.
#' @examples
#' x <- rnorm(20)
#' #--- All tests
#' (test <- symmetry_test(x))
#' #--- Choose any test
#' (Jk_test <- symmetry_test(x, statis = "Jk", Jk = 6))
#' #--- Choose several tests with multiple cut-offs
#' (MyTest1 <- symmetry_test(x, statis = c("Bk", "Jk"), Bk = c(5, 7), Jk = c(6, 8)))
#' (MyTest2 <- symmetry_test(x, statis = c("Mp", "Jk"), Jk = 6, Mp = c(10, 20, 25)))


symmetry_test <- function(x, statis = c("Bk", "Jk", "R", "Rs", "Mp", "Bkc", "Ck"), type = "k",
                          Bk = 5, Ck = 10, Jk = 6, Mp = c(10, 20, 25), median = 0) {
  
  if (type == "u") median <- median(x)
  
  n <- length(x)
  
  # --- Precálculos ---
  xc   <- x - median
  xabs <- abs(xc)
  si   <- (sign(xc) + 1) / 2
  ord  <- order(xabs)
  si   <- si[ord]
  indj <- seq_len(n)
  j    <- n - indj + 1
  Ij   <- c(1L, as.integer(diff(si) != 0))
  
  df <- data.frame(indj = indj, j = j, si = si, Ij = Ij)
  
  make_row <- function(stat, stat.value, p.value, n1 = NA_integer_, n0 = NA_integer_) {
    data.frame(stat       = stat,
               stat.value = round(stat.value, 1),
               p.value    = round(p.value, 4),
               n1         = n1,
               n0         = n0,
               stringsAsFactors = FALSE)
  }
  
  stats <- list()
  
  # --- Bk (múltiples cortes) ---
  if ("Bk" %in% statis) {
    for (k in Bk) {
      mask <- df$j <= k
      svBk <- sum(df$Ij[mask])
      stats[[paste0("Bk", k)]] <- make_row(paste0("Bk", k),
                                           svBk + 1,
                                           sum(stats::dbinom(0:svBk, k, 0.5)))
    }
  }
  

  # --- Ck (Bk condicional — versión robusta de Bk) ---

  if ("Ck" %in% statis) {
    for (k in Ck) {
      prim     <- n - k + 1
      mask     <- df$j < k
      svCk     <- sum(df$Ij[mask & df$indj >= 2])  # excluye Ij[1], igual que Bkc
      n1k      <- sum(df$si[mask]) + df$si[prim]   # signos positivos entre los k más cercanos
      n0k      <- k - n1k                     # signos negativos entre los k más cercanos
      stat.val <- svCk + 1                  # número de runs
      if (n1k == 0L || n0k == 0L) {
        pval_ck <- 1
      } else {
        pval_ck <- sum(druns(seq_len(stat.val) - 1, n1k, n0k))
      }
      stats[[paste0("Ck", k)]] <- make_row(paste0("Ck", k), stat.val, pval_ck, n1 = n1k, n0 = n0k)
    }
  }
  # --- Bkc: mismo estadistico que R, p-valor por distribucion condicional de runs ---

  if ("Bkc" %in% statis) {
    svBkc    <- sum(df$Ij[df$indj >= 2])   # identico a svR
    n1_bkc   <- sum(df$si)                  # signos positivos
    n0_bkc   <- n - n1_bkc                  # signos negativos
    stat.val <- svBkc + 1
    if (n1_bkc == 0L || n0_bkc == 0L) {
      pval_bkc <- 1
    } else {
      pval_bkc <- sum(druns(seq_len(stat.val) - 1, n1_bkc, n0_bkc)) # P(R* <= stat.val | n1, n0)
    }
    stats[["Bkc"]] <- make_row("Bkc", stat.val, pval_bkc, n1 = n1_bkc, n0 = n0_bkc)
  }
  
  # --- Jk (múltiples cortes) ---
  if ("Jk" %in% statis) {
    for (k in Jk) {
      mask <- df$j <= k
      svJk <- sum(df$indj[mask] * df$Ij[mask])
      EJk  <- k * (2*n - k + 1) / 4
      VJk  <- k * (6*n^2 + 6*n + 2*k^2 - 3*k - 6*n*k + 1) / 24
      stats[[paste0("Jk", k)]] <- make_row(paste0("Jk", k),
                                           svJk + 1,
                                           stats::pnorm((svJk - EJk) / sqrt(VJk)))
    }
  }
  
  # --- R ---
  if ("R" %in% statis) {
    svR <- sum(df$Ij[df$indj >= 2])
    stats[["R"]] <- make_row("R",
                             svR + 1,
                             sum(stats::dbinom(0:svR, n - 1, 0.5)))
  }
  
  # --- Rs ---
  if ("Rs" %in% statis) {
    svRs <- sum(df$Ij)
    n1   <- sum(df$si);  n0 <- n - n1
    ERs  <- 1 + (2*n0*n1) / n
    VRs  <- (2*n0*n1) * (2*n0*n1 - n) / (n^2 * (n - 1))
    stats[["Rs"]] <- make_row("Rs",
                              svRs,
                              stats::pnorm((svRs - ERs) / sqrt(VRs)))
  }
  
  # --- Mp ---
  if ("Mp" %in% statis) {
    for (p in Mp) {
      np   <- trunc(n * p / 100)
      mask <- df$indj >= np + 2
      svMp <- sum((df$indj[mask] - np) * df$Ij[mask])
      q    <- 1 - p/100
      EMp  <- (n*q - 1) * (n*q + 2) / 4
      VMp  <- (n*q - 1) * (2*n^2*q^2 + 5*n*q + 6) / 24
      stats[[paste0("M", p)]] <- make_row(paste0("M", p),
                                          svMp,
                                          stats::pnorm((svMp - EMp) / sqrt(VMp)))
    }
  }
  
  dplyr::bind_rows(stats)
}



