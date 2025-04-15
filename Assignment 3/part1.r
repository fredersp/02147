library(ggplot2)
library(patchwork)

library(httpgd)
hgd()


# Assignment 3: Part 1.1
set.seed(123)
rnorm(5)

phi_1 <- -0.6
phi_2 <- 0.5
n <- 200
lags <- 30

plot_simulations <- function(phi_1, phi_2) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()
  
  # Set up plotting layout for 5 plots in one figure
  par(mfrow = c(5, 1), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  
  # Function to plot each realization
  plotrealization <- function(x, i){
    plot(x, ylab = "X", main = paste("Simulation", i))
  }
  
  # Generate and plot 5 simulations
  for(i in 1:5){
    ts_data <- arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = 200)
    plotrealization(ts_data, i)
  }
}


plot_simulations_2 <- function(phi_1, phi_2) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()
  
  # Set up plotting layout for 5 plots in one figure
  par(mfrow = c(5, 1), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  
  # Function to plot each realization
  plotrealization <- function(x, i){
    plot(
      x,
      ylab = "X",
      xlab = "Time",
      main = paste("Simulation", i),
      cex.lab = 10,   # axis label size
      cex.axis = 10   # axis tick size
    )
  }
  
  # Generate and plot 5 simulations
  for(i in 1:5){
    ts_data <- arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = 200)
    plotrealization(ts_data, i)
  }
}

plot_simulations(phi_1, phi_2)
plot_simulations_2(phi_1, phi_2)




plot_simulations_gg <- function(phi_1, phi_2, n = 200) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()

  for (i in 1:5) {
    # Simulate AR(2) process
    ts_data <- arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = n)
    df <- data.frame(Time = 1:n, Value = ts_data)

    # Create plot with custom font size
    p <- ggplot(df, aes(x = Time, y = Value)) +
      geom_line() +
      labs(title = paste("Simulation", i), x = "Time", y = "X") +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
      )

    plot_list[[i]] <- p
  }

  # Stack all plots vertically
  final_plot <- wrap_plots(plot_list, ncol = 1)
  return(final_plot)
}

# Example usage:
plot_simulations_gg(phi_1, phi_2)
plot_simulations(phi_1, phi_2)




plot_acf_comparison <- function(phi_1, phi_2, n = 200, lags = 30) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()
  
  for (i in 1:5) {
    # Simulate AR(2)
    ts_data <- arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = n)
    
    # Empirical ACF
    emp_acf <- acf(ts_data, plot = FALSE, lag.max = lags)$acf
    
    # Theoretical ACF
    rho <- numeric(lags)
    rho[1] <- -phi_1 / (phi_2 + 1)
    rho[2] <- -phi_1 * rho[1] - phi_2
    for (k in 3:lags) {
      rho[k] <- -phi_1 * rho[k - 1] - phi_2 * rho[k - 2]
    }
    rho <- c(1, rho)
    
    # Combine into one data frame
    df <- data.frame(
      Lag = rep(0:lags, 2),
      ACF = c(emp_acf, rho),
      Type = rep(c("Empirical", "Theoretical"), each = lags + 1)
    )
    
    df$Lag_shifted <- df$Lag
    df$Lag_shifted[df$Type == "Empirical"] <- df$Lag_shifted[df$Type == "Empirical"] + 0.2  # Shift right by 0.2
    
    # Make plot
    p <- ggplot(df, aes(x = Lag_shifted, y = ACF, color = Type, fill = Type)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_segment(aes(xend = Lag_shifted, yend = 0), size = 1.2) +
      geom_point(size = 1.5) +
      scale_color_manual(values = c("Empirical" = "blue", "Theoretical" = "black")) +
      scale_fill_manual(values = c("Empirical" = "blue", "Theoretical" = "black")) +
      labs(title = paste("Empirical vs. Theoretical ACF — Realization", i),
           y = "ACF", x = "Lag") +
      scale_x_continuous(breaks = 0:lags) +
      theme_minimal()
    
    plot_list[[i]] <- p
  }
  
  # Combine into a 2x3 grid
  grid_plot <- (plot_list[[1]] | plot_list[[2]]) /
               (plot_list[[3]] | plot_list[[4]]) /
               (plot_list[[5]] + patchwork::plot_spacer())
  
  return(grid_plot)
}

library(ggplot2)
library(patchwork)

plot_acf_comparison <- function(phi_1, phi_2, n = 200, lags = 30) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()
  
  for (i in 1:5) {
    # Simulate AR(2)
    ts_data <- arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = n)
    
    # Empirical ACF
    emp_acf <- acf(ts_data, plot = FALSE, lag.max = lags)$acf
    
    # Theoretical ACF
    rho <- numeric(lags)
    rho[1] <- -phi_1 / (phi_2 + 1)
    rho[2] <- -phi_1 * rho[1] - phi_2
    for (k in 3:lags) {
      rho[k] <- -phi_1 * rho[k - 1] - phi_2 * rho[k - 2]
    }
    rho <- c(1, rho)
    
    # Combine into one data frame
    df <- data.frame(
      Lag = rep(0:lags, 2),
      ACF = c(emp_acf, rho),
      Type = rep(c("Empirical", "Theoretical"), each = lags + 1)
    )
    
    df$Lag_shifted <- df$Lag
    df$Lag_shifted[df$Type == "Empirical"] <- df$Lag_shifted[df$Type == "Empirical"] + 0.2  # Shift right
    
    # Make plot
    p <- ggplot(df, aes(x = Lag_shifted, y = ACF, color = Type, fill = Type)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_segment(aes(xend = Lag_shifted, yend = 0), size = 1.2) +
      geom_point(size = 1.5) +
      scale_color_manual(values = c("Empirical" = "blue", "Theoretical" = "black")) +
      scale_fill_manual(values = c("Empirical" = "blue", "Theoretical" = "black")) +
      labs(title = paste("Empirical vs. Theoretical ACF — Realization", i),
           y = "ACF", x = "Lag") +
      scale_x_continuous(breaks = seq(0, lags, by = 5)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12)
      )
    
    plot_list[[i]] <- p
  }
  
  # Combine vertically (1 column, 5 rows)
  grid_plot <- wrap_plots(plot_list, ncol = 1)
  return(grid_plot)
}

plot_acf_comparison(phi_1, phi_2, n = 200, lags = 30)


p <- plot_acf_comparison(phi_1, phi_2,n = 200, lags = 30)
ggsave("1.2 - acf_comparison_vertical.pdf", p, width = 8, height = 15)  # try height = 15 or higher

# 1.3

# Stationarty function
ar_roots <- function(phi_1, phi_2) {
  a <- 1
  b <- phi_1
  c <- phi_2

  discriminant <- b^2 - 4 * a * c

  # compute complex roots
  root1 <- (-b + sqrt(discriminant)) / (2 * a)
  root2 <- (-b - sqrt(discriminant)) / (2 * a)

  # calculate the modulus of the roots
  mod1 <- Mod(root1)
  mod2 <- Mod(root2)

  results <- list(mod1 = mod1, mod2 = mod2)
  return(results)

}


plot_simulations(-0.6, -0.3)
plot_acf_comparison(-0.6, -0.3, n = 200, lags = 30)
ar_roots(-0.6, -0.3)

# 1.4
plot_simulations(0.6, -0.3)
plot_acf_comparison(0.6, -0.3, n = 200, lags = 30)
ar_roots(0.6, -0.3)


# 1.5
plot_simulations(-0.7, -0.3)
plot_acf_comparison(-0.7, -0.3, n = 200, lags = 30)
ar_roots(-0.7, -0.3)


# 1.6
plot_simulations(-0.75, -0.3)
plot_acf_comparison(-0.75, -0.3, n = 200, lags = 30)
ar_roots(-0.75, -0.3)


# setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 3")