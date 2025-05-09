library(httpgd)
hgd() 


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 4")


# If needed, load your dataset
df <- read.csv("transformer_data.csv")

# Rename columns for clarity (optional, but helpful)
df <- df %>%
  rename(
    Yt = Y,
    Ta_t = Ta,
    Phi_s_t = S,
    Phi_I_t = I
  )


# Check structure
str(df)

# Reshape for plotting
df_long <- df %>%
  mutate(time = 1:n()) %>%
  pivot_longer(cols = c(Yt, Ta_t, Phi_s_t, Phi_I_t),
               names_to = "variable", values_to = "value")

# Plot all variables over time
ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(title = "Exploratory Time Series Plots",
       x = "Time (hours)",
       y = "Value") +
  theme_minimal()

