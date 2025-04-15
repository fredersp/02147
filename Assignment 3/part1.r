# Assignment 3 Part 1
library(ggplot2)
library(patchwork)

library(httpgd)
hgd()

 
# 1.1
set.seed(123)
rnorm(5)

phi_1 <- -0.6
phi_2 <- 0.5
n <- 200
lags <- 30


#### Simulation and plot of AR(2) process ####
plot_simulations_gg <- function(phi_1, phi_2, n = 200) {
  set.seed(123)
  rnorm(5)
  plot_list <- list()

  for (i in 1:5) {
    ts_data <- as.numeric(arima.sim(model = list(ar = c(-phi_1, -phi_2)), n = n))
    df <- data.frame(Time = 1:n, Value = ts_data)

    p <- ggplot(df, aes(x = Time, y = Value)) +
  geom_line(size = 0.4, color = "black") +
  labs(title = paste("Simulation", i), x = "Time", y = "X") +
  scale_y_continuous(breaks = seq(-3, 3, by = 2)) +  # Add ticks from -4 to 4
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 10),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    plot.margin = margin(5, 5, 5, 5)
  )

    plot_list[[i]] <- p
  }

  final_plot <- wrap_plots(plot_list, ncol = 1)
  return(final_plot)
}

# For phi_1 = -0.6 and phi_2 = 0.5:
plot_simulations_gg(phi_1, phi_2)



# 1.2

### Plot AFC comparision ####
plot_acf_comparison_gg <- function(phi_1, phi_2, n = 200, lags = 30) {
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

plot_acf_comparison_gg(phi_1, phi_2, n = 200, lags = 30)

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

# 1.2
plot_simulations_gg(-0.6, 0.5)
plot_acf_comparison_gg(-0.6, 0.5, n = 200, lags = 30)
#ar_roots(phi_1, phi_2)

# 1.3
plot_simulations_gg(-0.6, -0.3)
plot_acf_comparison_gg(-0.6, -0.3, n = 200, lags = 30)
ar_roots(-0.6, -0.3)

# 1.4
plot_simulations_gg(0.6, -0.3)
plot_acf_comparison_gg(0.6, -0.3, n = 200, lags = 30)
ar_roots(0.6, -0.3)


# 1.5
plot_simulations_gg(-0.7, -0.3)
plot_acf_comparison_gg(-0.7, -0.3, n = 200, lags = 30)
ar_roots(-0.7, -0.3)


# 1.6
plot_simulations_gg(-0.75, -0.3)
plot_acf_comparison_gg(-0.75, -0.3, n = 200, lags = 30)
ar_roots(-0.75, -0.3)


# setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 3")
