##' Plot an hnethist object
##'
##' Plots a heatmap or a BIC curve for an `hnethist` object.
##'
##' @param x an `hnethist` object from [hnethist()].
##' @param type One of `nethist` (default), `prob`, or `BIC`. When
##'   `type = "BIC"`, plots BIC values against the number of shapes `s`. The
##'   selected model is marked with a dashed vertical line, and the rightmost
##'   point (largest `s`, corresponding to the initial nethist fit) is labelled.
##'   `type = "bic"` is also accepted.
##' @param at A numeric vector of breakpoints for the color scale. Passed to
##'   [plot.nethist()]. See [plot.nethist()] for details.
##' @param ... Additional arguments passed to [plot()] when `type = "BIC"`, or
##'   to [plot.nethist()] otherwise.
##' @returns Called for its side effects (plotting). Returns `NULL` invisibly.
##' @seealso [plot.nethist()], [hnethist()]
##' @examples
##' \donttest{
##' set.seed(2022)
##' A <- igraph::as_adjacency_matrix(
##'   igraph::sample_gnp(100, 0.3), sparse = FALSE)
##' fit <- suppressMessages(hnethist(A))
##' plot(fit)
##' plot(fit, type = "BIC")
##' }
##' @importFrom graphics abline text
##' @exportS3Method
plot.hnethist <- function(x, type = "nethist", at = NULL, ...) {
  if (tolower(type) == "bic") {
    s_vals   <- sapply(x$details, function(d) d$s)
    bic_vals <- sapply(x$details, function(d) d$BIC)
    plot(s_vals, bic_vals, type = "b", xlab = "s", ylab = "BIC", ...)
    abline(v = x$s, lty = 2)
    text(max(s_vals), bic_vals[which.max(s_vals)], "initial", pos = 2)
    return(invisible(NULL))
  }
  plot.nethist(x, type = type, at = at, ...)
}

##' Network histogram plot
##'
##' Drawing [lattice::levelplot()] using a `nethist` object.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a
##'   permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`.
##' @param type One of `nethist` (default) or `prob`. `"pmat"` is a deprecated
##'   alias for `"prob"`.
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
##' @param at A numeric vector of breakpoints for the color scale. If `NULL`
##'   (default), breakpoints are computed automatically from the data range.
##'   Specify a fixed vector to compare plots from different fits on a common
##'   scale.
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
                         prob.col = "white", at = NULL, y = NA, ...) {
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

  at_vals    <- if (is.null(at)) seq(0, max(mat), length.out = 20) else at
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
##' @param type One of `nethist` (default) or `prob`. `"MNhist"` is a deprecated
##'   alias for `"nethist"`.
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
##' @param layout An integer vector `c(nrows, ncols)` specifying the panel
##'   grid for multi-layer plots, following the `mfrow` convention. If `NULL`
##'   (default), each layer is plotted as a separate figure.
##' @param layer_titles A character vector of length equal to the number of
##'   layers plotted, giving each panel's title. If `NULL` (default), titles
##'   default to `"Layer 1"`, `"Layer 2"`, etc.
##' @param at A numeric vector of breakpoints for the color scale. If `NULL`
##'   (default), breakpoints are computed automatically from the data range.
##'   Specify a fixed vector to compare plots from different fits on a common
##'   scale.
##' @param ... Additional arguments passed to [lattice::levelplot()].
##' @returns Called for its side effects (plotting). Returns `NULL` invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(IndianVil)
##' mnhist_Ind_vil <- multinethist(IndianVil)
##' plot(mnhist_Ind_vil)
##' plot(mnhist_Ind_vil, power = 0.5)
##' plot(mnhist_Ind_vil, layout = c(3,4))
##' plot(mnhist_Ind_vil, layer_titles = paste0("Network ", seq_along(mnhist_Ind_vil$rho_hat)))
##' }
##' @importFrom lattice levelplot panel.levelplot panel.text
##' @rdname plot.multinethist
##' @exportS3Method
plot.multinethist <- function(x, y = NA, type = "nethist",
                              idx_order = 1:max(x$cluster),
                              power = 0.25,
                              col.regions = function(n)
                                grDevices::hcl.colors(n, palette = "Reds 3", rev = TRUE),
                              colorkey = FALSE,
                              prob = FALSE, digits = 2,
                              prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
                              prob.col = "white",
                              layout = NULL,
                              layer_titles = NULL,
                              at = NULL, ...) {
  if (type == "MNhist") {
    warning("type = \"MNhist\" is deprecated. Use type = \"nethist\" instead.",
            call. = FALSE)
    type <- "nethist"
  }
  n_layers <- ifelse(length(dim(x$thetahat)) == 2, 1, dim(x$thetahat)[3])

  if (n_layers == 1) {
    return(invisible(plot.nethist(x, type, idx_order, power, col.regions,
                                  colorkey, prob, digits, prob.cex, prob.col,
                                  at = at, y = y, ...)))
  } else {
    return(invisible(.plot_multinethist_layers(x, type, idx_order, power,
                                               col.regions, colorkey, prob,
                                               digits, prob.cex, prob.col,
                                               layout = layout,
                                               layer_titles = layer_titles,
                                               at = at, ...)))
  }
}

.plot_multinethist_layers <- function(x, type = "nethist",
                                      idx_order = 1:max(x$cluster),
                                      power = 0.25,
                                      col.regions = function(n)
                                        grDevices::hcl.colors(n, palette = "viridis",
                                                              rev = TRUE),
                                      colorkey = FALSE,
                                      prob = FALSE, digits = 2,
                                      prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
                                      prob.col = "white",
                                      layout = NULL,
                                      layer_titles = NULL,
                                      at = NULL, ...) {
  k <- dim(x$thetahat)[1]
  if (!(type %in% c("nethist", "prob"))) {
    stop("type must be one of nethist or prob.")
  }
  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }

  # homogeneous case: only the first layer is shown
  n_loop <- ifelse(x$homogeneous & (type == "nethist"),
                   1L,
                   dim(x$thetahat)[3])

  if (!is.null(layout)) {
    if (!is.numeric(layout) || length(layout) != 2L || anyNA(layout) ||
        any(layout < 1L)) {
      stop("layout must be a length-2 integer vector of positive values.")
    }
    if (layout[1L] * layout[2L] < n_loop) {
      stop(paste0("layout c(", layout[1L], ", ", layout[2L], ") has only ",
                  layout[1L] * layout[2L], " cells but ", n_loop,
                  " layers to plot."))
    }
  }

  if (!is.null(layer_titles)) {
    if (!is.character(layer_titles) || length(layer_titles) != n_loop) {
      stop(paste0("layer_titles must be a character vector of length ", n_loop,
                  " (number of layers to plot)."))
    }
  }

  # compute global max for a consistent color scale across layers
  all_max <- vapply(seq_len(n_loop), function(l) {
    mat_l <- switch(type,
      "nethist" = (x$thetahat[idx_order, idx_order, l] / x$rho_hat[l])^power,
      "prob"   =  x$thetahat[idx_order, idx_order, l])
    max(mat_l)
  }, numeric(1L))
  at_vals <- if (is.null(at)) seq(0, max(all_max), length.out = 20) else at

  scales_arg <- list(
    x = list(at = seq_len(k), labels = as.character(idx_order)),
    y = list(at = seq_len(k), labels = as.character(idx_order))
  )

  extra <- list(...)

  if (!is.null(layout)) {
    nrows <- layout[1L]
    ncols <- layout[2L]
    if (!"aspect" %in% names(extra))
      extra[["aspect"]] <- "fill"
    if (!"par.settings" %in% names(extra))
      extra[["par.settings"]] <- list(
        layout.heights = list(top.padding = 0, bottom.padding = 0),
        layout.widths  = list(left.padding = 0, right.padding = 0)
      )
  }

  for (l in seq_len(n_loop)) {
    mat <- switch(type,
      "nethist" = (x$thetahat[idx_order, idx_order, l] / x$rho_hat[l])^power,
      "prob"   =  x$thetahat[idx_order, idx_order, l])
    rownames(mat) <- as.character(idx_order)
    colnames(mat) <- as.character(idx_order)

    default_title <- if (x$homogeneous && type == "nethist") NULL
                     else paste("Layer", l)
    base_args <- list(mat,
      col.regions = col.regions, at = at_vals, colorkey = colorkey,
      main = if (!is.null(layer_titles)) layer_titles[l] else default_title,
      ylim = c(k + 0.5, 0.5),
      xlab = "", ylab = "", scales = scales_arg)

    if (prob && type == "prob") {
      base_args[["panel"]] <- function(x, y, z, ...) {
        lattice::panel.levelplot(x, y, z, ...)
        lattice::panel.text(x, y, labels = round(z, digits),
                            col = prob.col, cex = prob.cex)
      }
    }
    p <- do.call(lattice::levelplot, c(base_args, extra))

    if (is.null(layout)) {
      print(p)
    } else {
      col_pos <- ((l - 1L) %% ncols) + 1L
      row_pos <- ((l - 1L) %/% ncols) + 1L
      print(p,
            split = c(col_pos, row_pos, ncols, nrows),
            more  = (l < n_loop))
    }
  }
  return(invisible(NULL))
}
