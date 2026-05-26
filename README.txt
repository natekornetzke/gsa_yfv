##This github repository contains code files for the manuscript 
##"Global Sensitivity Analysis for Stochastic Models of Infectious Diseases: Investigating the Re-Emergence of Yellow Fever Virus 
##in the Americas as a Case Study" by Nate Kornetzke and Helen J. Wearing
##Global sensitivity analysis code is contained in the R scripts
		-samples_for_simulation.R
		-sampling_and_sobol.R
		-modify_cpp_files.R
		-plotting_sobol_indices.R
##and the c++ files
		-run_A.cpp, run_B.cpp
		-run_C1.cpp, ..., run_C28.cpp
		-rates.cpp, rates.h
		-V.csv
##The two-at-a-time analysis code is contained in the zip files. All c++ was compiled using c++17 with gcc and openmp
	
