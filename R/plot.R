##' Network histogram plot
##'
##' Drawing [lattice::levelplot()] using a `nethist` object.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a
##'   permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`.
##' @param type One of `nethist` (default) or `prob`.
##' @param power A positive number for the power transform applied to the
##'   graphon estimate. Only used when `type = "nethist"`. Default is `0.25`.
##' @param col.regions A function taking an integer `n` and returning `n`
##'   colors. Default is viridis via [grDevices::hcl.colors()].
##' @param colorkey Logical. Whether to draw a color legend. Default `FALSE`.
##' @param prob Logical. Whether to print block probabilities on the plot.
##'   Default `FALSE`.
##' @param digits Integer. Number of decimal places for probabilities.
##' @param prob.cex Numeric. `cex` for probability labels.
##' @param prob.col Color for probability labels. Default `"white"`.
##' @param y A dummy variable for S3 dispatch. Never used.
##' @param ... Additional arguments passed to [lattice::levelplot()].
##' @returns Called for its side effects (plotting). Returns `NULL` invisibly.
##' @examples
##' \donttest{
##' set.seed(2022)
##' A <- igraph::sample_gnp(200, 0.05)
##' hist_A <- nethist(A)
##' plot(hist_A)
##' plot(hist_A, power = 0.5)
##' plot(hist_A, type = "prob", prob = TRUE)
##' }
##' @importFrom lattice levelplot panel.levelplot panel.text
##' @exportS3Method
plot.nethist <- function(x, type = "nethist",
                         idx_order = 1:max(x$cluster),
                         power = 0.25,
                         col.regions = function(n)
                           grDevices::hcl.colors(n, palette = "Reds 3", rev = TRUE),
                         colorkey = FALSE,
                         prob = FALSE, digits = 2,
                         prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
                         prob.col = "white", y = NA, ...) {
  k <- max(x$cluster)
  if (!(type %in% c("nethist", "prob", "pmat"))) {
    stop("type must be one of nethist or prob.")
  }
  if (type == "pmat") type <- "prob"

  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }

  mat <- switch(type,
    nethist = (x$thetahat[idx_order, idx_order] / x$rho_hat)^power,
    prob    =  x$thetahat[idx_order, idx_order])
  rownames(mat) <- as.character(idx_order)
  colnames(mat) <- as.character(idx_order)

  at_vals    <- seq(0, max(mat), length.out = 20)
  scales_arg <- list(
    x = list(at = seq_len(k), labels = as.character(idx_order)),
    y = list(at = seq_len(k), labels = as.character(idx_order))
  )

  if (prob && type == "prob") {
    p <- lattice::levelplot(mat,
          col.regions = col.regions, at = at_vals, colorkey = colorkey,
          ylim = c(k + 0.5, 0.5),
          xlab = "", ylab = "", scales = scales_arg,
          panel = function(x, y, z, ...) {
            lattice::panel.levelplot(x, y, z, ...)
            lattice::panel.text(x, y, labels = round(z, digits),
                                col = prob.col, cex = prob.cex)
          }, ...)
  } else {
    p <- lattice::levelplot(mat,
          col.regions = col.regions, at = at_vals, colorkey = colorkey,
          ylim = c(k + 0.5, 0.5),
          xlab = "", ylab = "", scales = scales_arg, ...)
  }
  print(p)
  return(invisible(NULL))
}

##' Network histogram plot
##'
##' Drawing [lattice::levelplot()] using a `multinethist` object.
##'
##' @param x a multinethist object from [multinethist()].
##' @param y A dummy variable for S3 dispatch. Never used.
##' @param idx_order A numeric vector for index label order, which must be a
##'   permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`.
##' @param type One of `MNhist` (default) or `prob`.
##' @param power A positive number for the power transform applied to the
##'   graphon estimate. Only used when `type = "MNhist"`. Default is `0.25`.
##' @param col.regions A function taking an integer `n` and returning `n`
##'   colors. Default is viridis via [grDevices::hcl.colors()].
##' @param colorkey Logical. Whether to draw a color legend. Default `FALSE`.
##' @param prob Logical. Whether to print block probabilities on the plot.
##'   Default `FALSE`.
##' @param digits Integer. Number of decimal places for probabilities.
##' @param prob.cex Numeric. `cex` for probability labels.
##' @param prob.col Color for probability labels. Default `"white"`.
##' @param ... Additional arguments passed to [lattice::levelplot()].
##' @returns Called for its side effects (plotting). Returns `NULL` invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(IndianVil)
##' mnhist_Ind_vil <- multinethist(IndianVil)
##' plot(mnhist_Ind_vil)
##' plot(mnhist_Ind_vil, power = 0.5)
##' }
##' @importFrom lattice levelplot panel.levelplot panel.text
##' @rdname plot.multinethist
##' @exportS3Method
plot.multinethist <- function(x, y = NA, type = "MNhist",
                              idx_order = 1:max(x$cluster),
                              power = 0.25,
                              col.regions = function(n)
                                grDevices::hcl.colors(n, palette = "Reds 3", rev = TRUE),
                              colorkey = FALSE,
                              prob = FALSE, digits = 2,
                              prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
                              prob.col = "white", ...) {
  n_layers <- ifelse(length(dim(x$thetahat)) == 2, 1, dim(x$thetahat)[3])

  if (n_layers == 1) {
    return(invisible(plot.nethist(x, type, idx_order, power, col.regions,
                                  colorkey, prob, digits, prob.cex, prob.col,
                                  y = y, ...)))
  } else {
    return(invisible(.plot_multinethist_layers(x, type, idx_order, power,
                                               col.regions, colorkey, prob,
                                               digits, prob.cex, prob.col, ...)))
  }
}

.plot_multinethist_layers <- function(x, type = "MNhist",
                                      idx_order = 1:max(x$cluster),
                                      power = 0.25,
                                      col.regions = function(n)
                                        grDevices::hcl.colors(n, palette = "viridis",
                                                              rev = TRUE),
                                      colorkey = FALSE,
                                      prob = FALSE, digits = 2,
                                      prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
                                      prob.col = "white", ...) {
  k <- dim(x$thetahat)[1]
  if (!(type %in% c("MNhist", "prob"))) {
    stop("type must be one of MNhist or prob.")
  }
  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }

  # homogeneous MNhist: only the first layer is shown
  n_loop <- ifelse(x$homogeneous & (type == "MNhist"),
                   1L,
                   dim(x$thetahat)[3])

  # compute global max for a consistent color scale across layers
  all_max <- vapply(seq_len(n_loop), function(l) {
    mat_l <- switch(type,
      "MNhist" = (x$thetahat[idx_order, idx_order, l] / x$rho_hat[l])^power,
      "prob"   =  x$thetahat[idx_order, idx_order, l])
    max(mat_l)
  }, numeric(1L))
  at_vals <- seq(0, max(all_max), length.out = 20)

  scales_arg <- list(
    x = list(at = seq_len(k), labels = as.character(idx_order)),
    y = list(at = seq_len(k), labels = as.character(idx_order))
  )

  for (l in seq_len(n_loop)) {
    mat <- switch(type,
      "MNhist" = (x$thetahat[idx_order, idx_order, l] / x$rho_hat[l])^power,
      "prob"   =  x$thetahat[idx_order, idx_order, l])
    rownames(mat) <- as.character(idx_order)
    colnames(mat) <- as.character(idx_order)

    if (prob && type == "prob") {
      print(lattice::levelplot(mat,
              col.regions = col.regions, at = at_vals, colorkey = colorkey,
              main = paste("Layer", l),
              ylim = c(k + 0.5, 0.5),
              xlab = "", ylab = "", scales = scales_arg,
              panel = function(x, y, z, ...) {
                lattice::panel.levelplot(x, y, z, ...)
                lattice::panel.text(x, y, labels = round(z, digits),
                                    col = prob.col, cex = prob.cex)
              }, ...))
    } else {
      print(lattice::levelplot(mat,
              col.regions = col.regions, at = at_vals, colorkey = colorkey,
              main = paste("Layer", l),
              ylim = c(k + 0.5, 0.5),
              xlab = "", ylab = "", scales = scales_arg, ...))
    }
  }
  return(invisible(NULL))
}
