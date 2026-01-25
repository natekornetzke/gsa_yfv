setwd("~/run_10")
source("~/run_10/sampling_and_sobol.R")

#number of samples and parameter space dimension
number_of_sample_realizations <- 3*10^5 + 5*10^4 # + 10^3
dim_of_problem <- 28

#intervals of parameters

#reference c++ model file to see which interval corresponds to which parameter
life_spans <- list(c(1.0/20.0, 1.0/12.0))
#life_spans <- list(c(1.0/20.0, 1.0/12.0),c(0.04, 0.08), c(0.05, 0.09), c(0.12, 0.16))
primate_birth <- list(c(0.1, 0.25))
mosquito_birth <- list(c(0.15, 0.45), c(0.11, 0.35), c(0.17, 0.47))
mosquito_biting <- list(c(21024, 46428), c(5383.75, 17976.25), c(11315, 30112.5))
human_to_mosquito_transmission <- list(c(0.08, 0.16), c(0.02, 0.2))
NHP_to_mosquito_transmission <- list(c(0.02, 0.2), c(0.2, 0.35))
mosquito_to_human_transmission <- list(c(0.25, 0.75), c(0.25, 0.75))
mosquito_to_NHP_transmission <- list(c(0.25, 0.75), c(0.25, 0.75))
primate_incubation <- list(c(0.008, 0.016), c(0.0014, 0.012))
mosquito_incubation <- list(c(0.019, 0.04), c(0.019, 0.04), c(0.019, 0.04))
recovery <- list(c(0.000391, 0.000913), c(0.000391, 0.000913))
primate_disease_death <- list(c(0.01, 0.62))

#human_travel <- list(c(0, 0.1), c(0, 0.1))
#primate_travel <- list(c(0,0.01),c(0,0.01))

human_travel <- list(c(0, 1), c(0, 1))
primate_travel <- list(c(0, 1),c(0, 1))

#previous simulations
#human_travel <- list(c(0, 0.001), c(0, 0.001))
#primate_travel <- list(c(0,0.001),c(0,0.001))
#best case
#human_travel <- list(c(0, 0.1), c(0, 0.1))
#primate_travel <- list(c(0,0.001),c(0,0.001))



#concatenate all the lists
parameter_intervals <- c(life_spans, primate_birth, mosquito_birth, mosquito_biting, human_to_mosquito_transmission, NHP_to_mosquito_transmission)
parameter_intervals <- c(parameter_intervals, mosquito_to_human_transmission, mosquito_to_NHP_transmission, primate_incubation)
#parameter_intervals <- c(parameter_intervals, mosquito_incubation, recovery, primate_disease_death, human_travel, primate_travel, carrying_capacities)
parameter_intervals <- c(parameter_intervals, mosquito_incubation, recovery, primate_disease_death, human_travel, primate_travel)
#double check to make sure we have the right number of parameter intervals
length(parameter_intervals)
if (length(parameter_intervals) != dim_of_problem)
{
  print("Error! Dimension mismatch.")
}

#generate sample matrices
sample_matrices <- generate_sample_matrices(dim_of_problem, number_of_sample_realizations, parameter_intervals)

#write sample matrices from list
A <- sample_matrices[[1]]
B <- sample_matrices[[2]]

write.table(A, file = "A.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(B, file = "B.csv", sep = ",", row.names = FALSE, col.names = FALSE)
for (i in 1:dim_of_problem) 
{
  file_name <- paste0("C", i, ".csv")
  write.table(sample_matrices[[i+2]], file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
}



