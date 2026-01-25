####Computing Sobol Indices for Stochastic Solvers
####Implementation based off of the mathematical framework presented in 
####Global Sensitivity Analysis in Monte Carlo Radiation Transport by Clements et al 2024 

compute_sobol_indices_deterministic <- function(y_A, y_B, y_C)
{
  S_i <- (mean(y_A*y_C) - mean(y_A)^2)/(mean(y_A^2) - mean(y_A)^2)
  T_i <- 1 - (mean(y_B*y_C) - mean(y_A)^2)/(mean(y_A^2) - mean(y_A)^2)
 
  indices_df <- data.frame("First Order" = S_i, "Total Order" = T_i)
  return(indices_df)
}


#function to compute QoI over samples, equation (7)
QoI_avg <- function(y_A, number_of_samples) 
{
  matrix_y <- matrix(y_A, nrow = number_of_samples)
  QoI_avg_result <- colMeans(matrix_y)
  return(QoI_avg_result)
}

#function to compute sampling estimator vector, equation (12)
solver_variances <- function(y_A, number_of_samples) 
{
  if (number_of_samples == 1)
  {
    solver_var_vec <- rep(0, length(y_A))
  }
  else
  {
    matrix_y <- matrix(y_A, nrow = number_of_samples)
    solver_var_vec <- apply(matrix_y, 2, var) 
  }
 
  return(solver_var_vec)
}
#function to compute polluted variance estimate, equation (14)
polluted_var <- function(y_A)
{
  polluted_S2 <- var(y_A)
  return(polluted_S2)
}

#function to compute parametric variance estimate, equation (15)
parametric_var <- function(solver_variances, polluted_S2, number_of_samples)
{
  number_of_sample_realizations <- length(solver_variances)
  S2 <- polluted_S2 - 1/(number_of_samples*number_of_sample_realizations)*sum(solver_variances)

  return(S2)
}

#function to compute first order and total order indices, using equations (16), (17), (18), (19)
compute_sobol_indices <- function(y_A, y_B, y_C, number_of_samples)
{
  #solver variances computed with respect to original model runs
  model_solver_variances <- solver_variances(y_A, number_of_samples)

  #other statistics computed with QoI average, equation (7)
  y_A_avg <- QoI_avg(y_A, number_of_samples)
  y_B_avg <- QoI_avg(y_B, number_of_samples)
  model_polluted_S2 <- polluted_var(y_A_avg)
  model_S2 <- parametric_var(model_solver_variances, model_polluted_S2, number_of_samples)
  
  y_C_i_avg <- QoI_avg(y_C, number_of_samples)
  S_i <- model_polluted_S2 - mean(y_A_avg^2 + y_C_i_avg^2 - (y_A_avg + y_C_i_avg)^2/2)
  S_i <- S_i/model_S2
  T_i <-  model_polluted_S2 - mean(y_B_avg^2 + y_C_i_avg^2 - (y_B_avg + y_C_i_avg)^2/2)
  T_i <- 1 - T_i/model_S2
  
  
  indices_df <- data.frame("First Order" = S_i, "Total Order" = T_i)
  
  
  return(indices_df)
}


####sampling tools, using the package "randtoolbox"
library(randtoolbox)

#function to scale intervals
scale_intervals <- function(A, intervals)
{
  dim <- length(intervals)
  for (i in 1:dim)
  {
    current_interval <-  intervals[[i]]
    a <-  current_interval[1]
    b <-  current_interval[2]
    A[,i] = a + (b-a)*A[,i]
  }
  return(A)
}

#function to generate matrices A, B, C_i using Sobol sequence, quasi-random numbers
generate_sample_matrices <- function(input_dimension, sample_size, intervals)
{
  #quasi-random sequence -- Sobol sequence
  sobol_seq <- sobol(n = sample_size, dim = 2*input_dimension)
  #grab A
  A <- sobol_seq[, 1:input_dimension]

  #grab B
  B <- sobol_seq[, (input_dimension + 1):(2 * input_dimension)]
  
  #scale matrices so we are sampling from the given intervals
  A <- scale_intervals(A, intervals)
  B <- scale_intervals(B, intervals)
 
  #generate C_i matrices
  # Initialize a list to store the matrices C_i
  sample_matrices <- vector("list", input_dimension)

  # Generate each matrix C_i
  for (i in 1:input_dimension) 
  {

    C_i <- B
    C_i[, i] <- A[, i]
    sample_matrices[[i]] <- C_i
  }
  
  #concatenate A, B to the front end of the samples_matrices vector
  sample_matrices <- c(list(A, B), sample_matrices)
  
  return(sample_matrices)
}

# ####test of sampling
# # n <- 10
# # d <- 3
# # interval <- c(-pi, pi)
# # intervals <-list(interval, interval, interval)
# # mats <- generate_sample_matrices(d, n, intervals)
# # mats
# 
# ####test case of sobol calculations with ishigami function
# ####note this is a test of a deterministic function
# ####see https://uqworld.org/t/ishigami-function/55
# ####or Sobol & Levitan 1999 for more details
# ishigami <- function(x1, x2, x3)
# {
#   ish <- sin(x1) + 7*sin(x2)^2 + 0.05*sin(x1)*x3^4
#   return(ish)
# }
# 
# 
# #sample size
# n <- 100000
# #dimension size
# d <- 3
# #qusi-random sequence -- Sobol sequence
# sobol_seq <- sobol(n = n, dim = 2*d)
# 
# #grab A
# A <- sobol_seq[, 1:d]
# #grab B
# B <- sobol_seq[, (d + 1):(2 * d)]
# 
# 
# interval <- c(-pi, pi)
# intervals = list(interval, interval, interval)
# A = scale_intervals(A, intervals)
# B = scale_intervals(B, intervals)
# 
# # Initialize a list to store the matrices C_i
# C_list <- vector("list", d)
# 
# # Generate each matrix C_i
# for (i in 1:d) {
#   C_i <- B
#   C_i[, i] <- A[, i]
#   file_name <- paste0("C", i, ".csv")
#   #write.csv(C_i, file = file_name, row.names = FALSE)
#   C_list[[i]] <- C_i
# }
# 
# C1 <- C_list[[1]]
# C2 <- C_list[[2]]
# C3 <- C_list[[3]]
# 
# ish_A <- ishigami(A[,1], A[,2], A[,3])
# ish_B <- ishigami(B[,1], B[,2], B[,3])
# ish_C1 <- ishigami(C1[,1], C1[,2], C1[,3])
# ish_C2 <- ishigami(C2[,1], C2[,2], C1[,3])
# ish_C3 <- ishigami(C3[,1], C3[,2], C3[,3])
# y_Cs <- vector(mode = "list", 3)
# y_Cs[[1]] <- ish_C1
# y_Cs[[2]] <- ish_C2
# y_Cs[[3]] <- ish_C3
# 
# indices <- c(1, 2, 3)
# 
# 
# sobol_indices <- compute_sobol_indices(ish_A, ish_B, y_Cs, indices, 1)
# sobol_indices


