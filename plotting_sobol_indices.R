setwd("~/Desktop/runs_3_9_25/march_2_2025/")
library(ggplot2)
library(tidyr)
library(dplyr)



# enum ParamIndices {
#   mu_p, vu_p,
#   Theta_u, Theta_p, Theta_s,
#   b_u, b_p, b_s,
#   beta_hm_u, beta_hm_p,
#   beta_pm_p, beta_pm_s,
#   beta_mh_u, beta_mh_p,
#   beta_mp_p, beta_mp_s,
#   epsilon_h, epsilon_p,
#   epsilon_u_m, epsilon_p_m, epsilon_s_m,
#   rho_h, rho_p,
#   eta_p,
#   tau_u, tau_p,
#   varphi_p, varphi_s,
# };
x_labels <- c(expression(mu[p]), expression(vu[p]), expression(Theta[u]), expression(Theta[p]), expression(Theta[s]))
x_labels <- c(x_labels, expression(b[u]), expression(b[p]), expression(b[s]), expression(beta^(hm)[u]), expression(beta^(hm)[p]))
x_labels <- c(x_labels, expression(beta^{pm}[p]), expression(beta^(pm)[s]), expression(beta^(mh)[u]), expression(beta^(mh)[p]))
x_labels <- c(x_labels, expression(beta^(mp)[p]), expression(beta^(mp)[s]), expression(epsilon[h]), expression(epsilon[p]))
x_labels <- c(x_labels, expression(epsilon^u[m]), expression(epsilon^p[m]), expression(epsilon^s[m]), expression(rho[h]), expression(rho[p]))
x_labels <- c(x_labels, expression(eta[p]), expression(tau[u]), expression(tau[p]), expression(varphi[p]), expression(varphi[s]))
length(x_labels)

df_SNHP <- read.csv("SNHP_indices.csv")
df_SNHP$First.Order[df_SNHP$First.Order < 0] <- 0 
df_SNHP$First.Order[df_SNHP$First.Order < 0.01] <- NA 
df_SNHP$Total.Order[df_SNHP$Total.Order < 0] <- 0 
df_SNHP$Total.Order[df_SNHP$Total.Order < 0.01] <- NA 
df_PNHP <- read.csv("PNHP_indices.csv")
df_PNHP$First.Order[df_PNHP$First.Order < 0] <- 0
df_PNHP$First.Order[df_PNHP$First.Order < 0.01] <- NA
df_PNHP$Total.Order[df_PNHP$Total.Order < 0] <- 0
df_PNHP$Total.Order[df_PNHP$Total.Order < 0.01] <- NA
df_PH <- read.csv("PH_indices.csv")
df_PH$First.Order[df_PH$First.Order < 0] <- 0
df_PH$First.Order[df_PH$First.Order < 0.01] <- NA
df_PH$Total.Order[df_PH$Total.Order < 0] <- 0
df_PH$Total.Order[df_PH$Total.Order < 0.01] <- NA
df_UH <- read.csv("UH_indices.csv")
df_UH$First.Order[df_UH$First.Order < 0] <- 0
df_UH$First.Order[df_UH$First.Order < 0.01] <- NA
df_PH$Total.Order[df_PH$Total.Order < 0] <- 0
df_PH$Total.Order[df_PH$Total.Order < 0.01] <- NA

first_order_indices <- data.frame(SNHP = df_SNHP$First.Order, PNHP = df_PNHP$First.Order, PH = df_PH$First.Order, UH = df_UH$First.Order)
total_order_indices <- data.frame(SNHP = df_SNHP$Total.Order, PNHP = df_PNHP$Total.Order, PH = df_PH$Total.Order, UH = df_UH$Total.Order)


#use tidyr to graph
first_order_indices_long <- first_order_indices %>%
  pivot_longer(
    cols = everything(),    # Include all columns to pivot
    names_to = "QOI",     # New column for original column names
    values_to = "Value"      # New column for values
  )
first_order_indices_long <- first_order_indices_long %>% mutate(Param = rep(1:28, each = 4))

total_order_indices_long <- total_order_indices %>%
  pivot_longer(
    cols = everything(),    # Include all columns to pivot
    names_to = "QOI",     # New column for original column names
    values_to = "Value"      # New column for values
  )
total_order_indices_long <- total_order_indices_long %>% mutate(Param = rep(1:28, each = 4))



#first order indices for each QoI
p1 <- ggplot(data = first_order_indices_long, mapping = aes(x = as.factor(Param), y = QOI, fill = Value))
p1 <- p1 + geom_tile(color = "black") + scale_fill_gradient(low = "skyblue", high = "red", na.value = "white")
p1 <- p1 + theme_minimal() +   theme(
  panel.background = element_blank(),       # Removes the gray background
  panel.grid = element_blank(),             # Removes all grid lines
  axis.line = element_line(color = "black") # Adds black lines for the x and y axes
)
p1 <- p1 + scale_x_discrete(labels = x_labels, expand = c(0,0)) + scale_y_discrete(expand = c(0,0))
p1 <- p1 + labs(x = "Parameters", y = "Quantities of Interest", fill = "First Order Sobol Index")
p1

#total order indices for each QoI
p2 <- ggplot(data = total_order_indices_long, mapping = aes(x = as.factor(Param), y = QOI, fill = Value))
p2 <- p2 + geom_tile(color = "black") + scale_fill_gradient(low = "skyblue", high = "red", na.value = "white")
p2 <- p2 + theme_minimal() +   theme(
  panel.background = element_blank(),       # Removes the gray background
  panel.grid = element_blank(),             # Removes all grid lines
  axis.line = element_line(color = "black") # Adds black lines for the x and y axes
)
p2 <- p2 + scale_x_discrete(labels = x_labels, expand = c(0,0)) + scale_y_discrete(expand = c(0,0))
p2 <- p2 + labs(x = "Parameters", y = "Quantities of Interest", fill = "Total Order Sobol Index")
p2