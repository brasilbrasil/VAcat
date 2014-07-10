#' Convert a continuous vector to discrete values
#'
#' Take a vector of numeric values and return a vector of categorical values
#' based on either quantiles or absolute values. If there are multiple
#' duplicates at the extremes so that a cutoff is "overrun" (e.g., 50% zeros
#' with a lower cutpoint of 0.33) the cutpoint is set to be just outside the
#' minimum or maximum value.
#' @param x vector of numeric values
#' @param q vector of (2) cutpoints for the quantile function. Defaults to
#' c(0.33, 0.66)
#' @param lv a character vector of (3) levels to name the categories. Defaults
#' c("Low", "Medium", "High")
#' @param absolute if TRUE, then the cutpoints in q are taken as
#' given. If FALSE they are used as as cutpoints for quantiles of vector x.
#' Defaults to TRUE.
#' @param eps the small amount to adjust the cutpoints in the case of
#' overrun as described above. Defaults to 0.0001
#' @export
q3levels <- function(x, q=c(0.33, 0.66), lv=c("Low", "Medium", "High"),
                     absolute=TRUE, eps=0.0001) {
    if (!absolute) {
        y <- quantile(x, q, na.rm=TRUE)
        if (y[1] == min(x, na.rm=TRUE)) {
            y[1] <- y[1] + eps
        }
        if (y[2] == max(x, na.rm=TRUE)) {
            y[2] <- y[2] - eps
        }
    }
    quant_data <- rep(lv[2], length(x))
    quant_data[is.na(x)] <- "none"
    quant_data[x < y[1]] <- lv[1]
    quant_data[x > y[2]] <- lv[3]
    return(list(quant_data, thresholds=y))
}

#' Convert a continuous vector to discrete values
#'
#' Take a vector of numeric values and return a vector of categorical values
#' based on either quantiles or absolute values. If there are multiple
#' duplicates at the extremes so that a cutoff is "overrun" (e.g., 50% zeros
#' with a lower cutpoint of 0.33) the cutpoint is set to be just outside the
#' minimum or maximum value.
#' @param x vector of numeric values
#' @param q a single cutpoint for the quantile function. Defaults to 0.5.
#' @param lv a character vector of (2) levels to name the categories. Defaults
#' c("Low", "High")
#' @param absolute if TRUE, then the cutpoint given in q is taken as
#' given. If FALSE they it is used as the cutpoint for vector x. Defaults
#' to TRUE
#' @param eps the small amount to adjust the cutpoint in the case of
#' overrun as described above. Defaults to 0.0001
#' @export
q2levels <- function(x, q=0.5, lv=c("Low", "High"),
                     absolute=TRUE, eps=0.0001) {
    if (!absolute) {
        y <- quantile(x, q, na.rm=TRUE)
        if (y == min(x, na.rm=TRUE)) {
            y <- y + eps
        } else if (y == max(x, na.rm=TRUE)) {
            y <- y - eps
        }
    }
  quant_data <- rep(lv[1], length(x))
  quant_data[is.na(x)] <- "none"
  quant_data[x > y] <- lv[2]
  return(list(quant_data, thresholds=y))
}

#' Convert a continuous vector to discrete values
#'
#' Take a vector of numeric values and return a vector of categorical values
#' based on either quantiles or absolute values. If there are multiple
#' duplicates at the extremes so that a cutoff is "overrun" (e.g., 50% zeros
#' with a lower cutpoint of 0.33) the cutpoint is set to be just outside the
#' minimum or maximum value.
#' @param x vector of numeric values
#' @param q vector of (4) cutpoints for the quantile function. Defaults to 20%
#' in each category
#' @param lv a character vector of (5) levels to name the categories. Defaults
#' range from "Very_small" to "Very_large".
#' @param absolute if TRUE, then the cutpoints in q are taken as
#' given. If FALSE they are used as as cutpoints for quantiles of vector x.
#' Defaults to TRUE.
#' @param eps the small amount to adjust the cutpoints in the case of
#' overrun as described above. Defaults to 0.0001
#' @export
q5levels <- function(x, q=c(0.2, 0.4, 0.6, 0.8),
                     lv=c("Very_small", "Small", "Medium", "Large", "Very_large"),
                     absolute=TRUE, eps=0.0001) {
    if (!absolute) {
        y <- quantile(x, q, na.rm=TRUE)
        if (y[1] == min(x, na.rm=TRUE)) {
            y[1] <- y[1] + eps
        }
        if (y[4] == max(x, na.rm=TRUE)) {
            y[4] <- y[4] - eps
        }
    }
    quant_data <- rep(lv[3], length(x))
    quant_data[is.na(x)] <- "none"
    quant_data[x < y[2]] <- lv[2]
    quant_data[x < y[1]] <- lv[1]
    quant_data[x > y[3]] <- lv[4]
    quant_data[x > y[4]] <- lv[5]
    return(list(quant_data, thresholds=y))
}
