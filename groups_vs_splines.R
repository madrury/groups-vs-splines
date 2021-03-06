library(ggplot2)
library(splines)

test_cuts_vs_splines <- function(signal, N, noise,
                                 range=c(0, 1),
                                 max_parameters=50,
                                 seed=154) {

  if(max_parameters < 8) {
    stop("Please pass max_parameters >= 8, otherwise the plots look kinda bad.")
  }

  out_obj <- list()

  set.seed(seed)

  #---- Make training and testing data sets.
  x_train <- runif(N, range[1], range[2])
  x_test <- runif(N, range[1], range[2])

  y_train <- signal(x_train) + rnorm(N, 0, noise)
  y_test <- signal(x_test) + rnorm(N, 0, noise)

  out_obj$training_data <- data.frame(x = x_train, y = y_train)
  out_obj$testing_data <- data.frame(x = x_test, y = y_test)

  #---- Draw a simple line plot of the true signal
  df_signal_plot <- data.frame(
    x = seq(range[1], range[2], length.out = 100)
  )
  df_signal_plot$y <- signal(df_signal_plot$x)
  out_obj$signal_plot <- ggplot(data = df_signal_plot) +
    geom_line(aes(x = x, y = y)) +
    labs(title = "True Signal")

  #---- A plot of the training and testing data
  df_data_plot <- data.frame(
    x = c(x_train, x_test),
    y = c(y_train, y_test),
    id = c(rep("train", N), rep("test", N))
  )
  out_obj$data_plot <- ggplot(data = df_data_plot) +
    geom_point(aes(x=x, y=y)) +
    facet_wrap(~ id) +
    labs(title = "Training and Testing Data")

  #--- A plot of the signal and data together
  out_obj$signal_and_data_plot <- ggplot() +
    geom_line(aes(x = x, y = y), color = 'blue', size = 2, alpha = 0.6,
              data = df_signal_plot) +
    geom_point(aes(x = x, y = y), data = df_data_plot) +
    facet_wrap(~ id) +
    labs(title = "Training and Testing Data with True Signal")

  #----- lm with various groupings -------------
  models_with_groupings <- list()
  train_errors_cuts <- rep(NULL, length(models_with_groupings))
  test_errors_cuts <- rep(NULL, length(models_with_groupings))

  for (n_groups in 3:max_parameters) {
    cut_points <- seq(range[1], range[2], length.out = n_groups + 1)
    x_train_factor <- cut(x_train, cut_points)
    factor_train_data <- data.frame(x = x_train_factor, y = y_train)
    models_with_groupings[[n_groups]] <- lm(y ~ x, data = factor_train_data)

    # Training error rate
    train_preds <- predict(models_with_groupings[[n_groups]], factor_train_data)
    soses <- (1/N) * sum( (y_train - train_preds)**2)
    train_errors_cuts[n_groups - 2] <- soses

    # Testing error rate
    x_test_factor <- cut(x_test, cut_points)
    factor_test_data <- data.frame(x = x_test_factor, y = y_test)
    test_preds <- predict(models_with_groupings[[n_groups]], factor_test_data)
    soses <- (1/N) * sum( (y_test - test_preds)**2)
    test_errors_cuts[n_groups - 2] <- soses
  }
  out_obj$models_with_groupings <- models_with_groupings

  error_df_cuts <- data.frame(
    x = rep(3:max_parameters, 2),
    e = c(train_errors_cuts, test_errors_cuts),
    id = c(rep("train", length(train_errors_cuts)),
           rep("test", length(test_errors_cuts))),
    type = "cuts"
  )
  out_obj$errors_cuts_plot <- ggplot(data = error_df_cuts) +
    geom_line(aes(x = x, y = e)) +
    facet_wrap(~ id) +
    labs(title = "Error Rates with Grouping Transformations",
         x = ("Number of Estimated Parameters"),
         y = ("Average Squared Error"))

  #----- lm with natural splines -------------
  models_with_splines <- list()
  train_errors_splines <- rep(NULL, length(models_with_groupings))
  test_errors_splines <- rep(NULL, length(models_with_groupings))

  for (deg_freedom in 3:max_parameters) {
    knots <- seq(range[1], range[2], length.out = deg_freedom + 1)[2:deg_freedom]

    train_data <- data.frame(x = x_train, y = y_train)
    models_with_splines[[deg_freedom]] <- lm(y ~ ns(x, knots=knots), data = train_data)

    # Training error rate
    train_preds <- predict(models_with_splines[[deg_freedom]], train_data)
    soses <- (1/N) * sum( (y_train - train_preds)**2)
    train_errors_splines[deg_freedom - 2] <- soses

    # Testing error rate
    test_data <- data.frame(x = x_test, y = y_test)
    test_preds <- predict(models_with_splines[[deg_freedom]], test_data)
    soses <- (1/N) * sum( (y_test - test_preds)**2)
    test_errors_splines[deg_freedom - 2] <- soses
  }
  out_obj$models_with_splines <- models_with_splines

  error_df_splines <- data.frame(
    x = rep(3:max_parameters, 2),
    e = c(train_errors_splines, test_errors_splines),
    id = c(rep("train", length(train_errors_splines)),
           rep("test", length(test_errors_splines))),
    type = "splines"
  )
  out_obj$errors_splines_plot <- ggplot(data = error_df_splines) +
    geom_line(aes(x = x, y = e)) +
    facet_wrap(~ id) +
    labs(title = "Error Rates with Natural Cubic Spline Transformations",
         x = ("Number of Estimated Parameters"),
         y = ("Average Squared Error"))


  error_df <- rbind(error_df_cuts, error_df_splines)
  out_obj$error_df <- error_df

  # The training error for the first cut model is always an outlier, and
  # messes up the y range of the plots.
  y_lower_bound <- min(c(train_errors_cuts, train_errors_splines))
  y_upper_bound = train_errors_cuts[2]
  out_obj$errors_comparison_plot <- ggplot(data = error_df) +
    geom_line(aes(x = x, y = e)) +
    facet_wrap(~ id*type) +
    scale_y_continuous(limits = c(y_lower_bound, y_upper_bound)) +
    labs(
      title = ("Binning vs. Natural Splines"),
      x = ("Number of Estimated Parameters"),
      y = ("Average Squared Error"))

  out_obj
}


#---- A very non linear function
true_signal_weird <- function(x) {
  x*x*x*(x-1) + 2*(1/(1+exp(-.5*(x-.5)))) - 3.5*(x > .2)*(x < .5)*(x - .2)*(x - .5)
}
obj <- test_cuts_vs_splines(true_signal_weird, 250, 0.075)


# Plot the data and predicted values from the grouping model
plot_data_and_predictions_cuts <- function(obj, n_groups) {
  M <- obj$models_with_groupings[[n_groups]]

  plot_df <- data.frame(x_cont=seq(0, 1, 0.005))
  plot_df$x <- cut(plot_df$x_cont, seq(0, 1, length.out = n_groups + 1))
  plot_df$y <- predict(M, newdata=plot_df)

  ggplot() +
    geom_line(aes(x=x_cont, y=y), data=plot_df, color="blue", size=2, alpha=0.5) +
    geom_point(aes(x=x, y=y), data=obj$testing_data)
}
plot_data_and_predictions_cuts(obj, 5)

# Plot the data and predicted values from the spline model
plot_data_and_predictions_splines <- function(obj, n_knots) {

  M <- obj$models_with_splines[[n_knots]]

  plot_df <- data.frame(x=seq(0, 1, 0.005))
  plot_df$y <- predict(M, newdata=plot_df)

  ggplot() +
    geom_line(aes(x=x, y=y), data=plot_df, color="blue", size=2, alpha=0.5) +
    geom_point(aes(x=x, y=y), data=obj$testing_data)
}
plot_data_and_predictions_splines(obj, 6)

#--- Let's run the experiment a ton of times to make sure of our conclusions
replicate_error_comparison <- function(signal, n_replicates, N=500, noise=0.05) {
  replicates <- list()
  for(i in 1:n_replicates) {
    replicates[[i]] <- test_cuts_vs_splines(signal, N, noise, seed=100*i)
  }

  K <- nrow(replicates[[1]]$error_df)
  errors <- do.call(rbind, lapply(replicates, `[[`, "error_df"))
  errors$model_id <- rep(1:n_replicates, each=K)

  # Summarize and make the average error curve
  errors_summary <- errors %>%
    group_by_(.dots=c("x", "id", "type")) %>%
    summarize(e = mean(e))

  ggplot() +
    geom_line(aes(x=x, y=e, group=model_id), alpha=0.25, data=errors) +
    geom_line(aes(x=x, y=e), data=errors_summary, size=2) +
    facet_wrap(~ id*type)
}
replicate_error_comparison(true_signal_weird, 50)

#--- A sinusoidal signal function

true_signal_sin <- function(x) {
  x + 1.5*sin(3*2*pi*x)
}
obj <- test_cuts_vs_splines(true_signal_sin, 250, 1)
replicate_error_comparison(true_signal_sin, 50, 500, 1)
plot_data_and_predictions_cuts(obj, 25)
plot_data_and_predictions_splines(obj, 10)

true_signal_line <- function(x) {x}
obj <- test_cuts_vs_splines(true_signal_line, 250, .2)
replicate_error_comparison(true_signal_line, 50, 500, 0.2)
plot_data_and_predictions_cuts(obj, 12)
plot_data_and_predictions_splines(obj, 20)
