
#' Plot predictions of a gam model
#'
#' Prediction plot of a gam model, including confidence intervals for the mean.
#' Also, the ggplot2 object and a table of the predictions is silently returned as a list.
#' @param model gam model
#' @param nLetters.WD a numeric value specifying the constant \code{nLetters.WD}
#' value that should be used for all observations in newdat (original values will be overridden).
#' Alternatively the value \code{"use_mean"} can be passed. In this case
#' the mean observed \code{nLetters.WD} value is used for each AOI.condition.
#' @param groups_vec character vector of categories of the grouping variable (usually \code{"AOI.condition"})
#' to be plotted.
#' @param group_means Optional data frame with columns \code{"AOI.condition"} and
#' \code{nLetters.WD} to pass group specific means for predictions.
#' Only works in combination with \code{nLetters.WD = "use_mean"}.
#' @param model_name name used in the plot title
#' @param grouping_var name of the grouping variable. Defaults to \code{"AOI.condition"}.
#' This argument can e.g. be used when also language groups are included in the grouping
#' and the grouping variable is \code{"AOI.condition.group"}.
#' @param nonnegative_predictions If \code{TRUE} (default), negative predictions
#' of the mean value are set to zero.
#' @param xlab_size size of text at x axis ticks
#' @param point_size size of the points in the plot
#' @param ylim limits for y axis
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @export
#' @return list containing the ggplot2 object and a table of predictions
plot_predictions <- function(model, nLetters.WD, groups_vec, group_means = NULL,
                             model_name = "MODEL_NAME", grouping_var = "AOI.condition",
                             nonnegative_predictions = TRUE,
                             xlab_size = 6, point_size = 2, ylim = waiver()) {
  if (missing(nLetters.WD) | (!is.numeric(nLetters.WD) && nLetters.WD != "use_mean"))
    stop("Please specify the argument 'nLetters.WD' according to ?plot_predictions")

  response <- as.character(model$formula[2])
  response_short <- substring(response, 1, 3)

  # give the grouping variable a standard name
  model_dat <- model$model
  colnames(model_dat)[colnames(model_dat) == grouping_var] <- "group"

  # calculate mean nLetters.WD per group
  mean_dat <- model_dat %>% filter(group %in% groups_vec) %>%
    group_by(group) %>% summarize(nLetters.WD = mean(nLetters.WD)) %>%
    mutate(group = as.character(group))

  ### create the new groupings dataset for the predictions
  if (nLetters.WD == "use_mean") {
    if (!is.null(group_means)) {
      newdat <- group_means %>% filter(group %in% groups_vec)
    } else {
      newdat <- mean_dat
    }
    attr(newdat, "terms") <- NULL
  } else if (is.numeric(nLetters.WD)) {
    newdat <- data.frame("group" = unique(model_dat$group),
                         "nLetters.WD"   = nLetters.WD) %>%
      filter(group %in% groups_vec)
  }
  # include random participants and topics (as these are excluded anyway in the predictions)
  newdat$Participant <- sample(model_dat$Participant, 1)
  newdat$Topic <- sample(model_dat$Topic, 1)

  # calculate predictions
  colnames(newdat)[colnames(newdat) == "group"] <- grouping_var
  p <- predict.gam(model, newdat, se.fit = TRUE, exclude = c("s(Participant)","s(Topic)"))
  colnames(newdat)[colnames(newdat) == grouping_var] <- "group"

  ### create a table with the predictions and confidence intervals
  pred_dat <- newdat %>%
    select(group, nLetters.WD) %>%
    dplyr::rename(nLetters.WD_fix = nLetters.WD) %>%
    mutate(group = as.character(group)) %>%
    left_join(mean_dat, by = "group") %>%
    dplyr::rename(nLetters.WD_obs = nLetters.WD) %>%
    mutate(group = factor(as.character(group), levels = groups_vec)) %>%
    mutate(pred = p$fit,
           pred_StdErr = p$se.fit,
           ci_lower = p$fit - qnorm(0.975) * p$se.fit,
           ci_upper = p$fit + qnorm(0.975) * p$se.fit)

  ### add the model estimates
  est <- as.data.frame(summary(model)$p.table)
  # if not all groups appear in row.names(est), the missing one is the intercept.
  # before applying grepl, remove '+' signs as these lead to problems
  groups_vec_check <- gsub("\\+", "", groups_vec)
  rn_est_check <- gsub("\\+", "", row.names(est))
  check <- sapply(groups_vec_check, function(x) any(grepl(x, rn_est_check)))
  if (!all(check))
    row.names(est)[1] <- paste0(grouping_var, groups_vec[!check])
  # add estimates and standard errors to our table
  est$group <- factor(substr(row.names(est), nchar(grouping_var) + 1, 50), levels = groups_vec)
  est <- est %>% filter(!is.na(est$group))
  pred_dat <- pred_dat %>%
    left_join(est[,c("Estimate","Std. Error","group")], by = "group")
  colnames(pred_dat)[colnames(pred_dat) == "Std. Error"] <- "StdErr"

  # reshape data for plotting
  plot_dat <- pred_dat %>%
    select(group, pred, ci_lower, ci_upper) %>%
    gather(key = type, value = value, -group)
  # specific preparations for drawing the horizontal lines at the borders of the confidence intervals
  plot_dat2 <- bind_rows(plot_dat, plot_dat) %>%
    mutate(border_id = paste0(group, type)) %>%
    mutate(x_num = as.numeric(group))
  plot_dat2$x_num[1:(nrow(plot_dat2)/2)] <- plot_dat2$x_num[1:(nrow(plot_dat2)/2)] - 0.2
  plot_dat2$x_num[(nrow(plot_dat2)/2+1):nrow(plot_dat2)] <- plot_dat2$x_num[(nrow(plot_dat2)/2+1):nrow(plot_dat2)] + 0.2
  plot_dat2 <- plot_dat2 %>% filter(type != "pred")

  ### set negative values to zero
  if (nonnegative_predictions) {
    pred_dat$Estimate[pred_dat$pred < 0] <- 0
    pred_dat$pred[pred_dat$pred < 0]     <- 0
    plot_dat$value[plot_dat$value < 0 & plot_dat$type == "pred"] <- 0
  }
  
  ### further plot preparations
  if (!identical(ylim, waiver())) {
    # extend ylim if the confidence interval borders are wider than the ylim borders
    if (ylim[1] > min(plot_dat$value))
      ylim[1] <- min(plot_dat$value)
    if (ylim[2] < max(plot_dat$value))
      ylim[2] <- max(plot_dat$value)
  }

  ### plot
  gg <- ggplot() +
    geom_line(data = plot_dat %>% filter(type != "pred"), aes(x = group, y = value, group = group)) +
    geom_line(data = plot_dat2, aes(x = x_num, y = value, group = border_id)) +
    geom_point(data = plot_dat %>% filter(type == "pred"), aes(x = group, y = value), size = point_size) +
    xlab(grouping_var) + ylab(paste(response_short, "per word")) +
    ggtitle(paste(model_name, "-", response_short)) +
    scale_y_continuous(limits = ylim) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = xlab_size),
          plot.title = element_text(hjust = 0.5))
  print(gg)

  ### final preparations of the table
  pred_dat <- pred_dat %>%
    select(group, Estimate, StdErr, nLetters.WD_obs, nLetters.WD_fix, pred, pred_StdErr)
  colnames(pred_dat)[colnames(pred_dat) == "pred"] <- paste0(response_short, ".Pred")
  colnames(pred_dat)[colnames(pred_dat) == "pred_StdErr"] <- paste0(response_short, ".Pred.StdErr")
  pred_dat <- as.data.frame(pred_dat)
  row.names(pred_dat) <- pred_dat$group
  pred_dat$group <- NULL
  pred_dat <- as.matrix(pred_dat)
  pred_dat <- pred_dat[match(groups_vec, row.names(pred_dat)),]
  pred_dat <- round(pred_dat, 2)

  # invisibly return the table
  invisible(list("plot" = gg, "table" = pred_dat))
}
