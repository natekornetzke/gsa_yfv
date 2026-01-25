setwd("~/Desktop/binomial_tau/")

# Read the file
lines <- readLines("rates.txt")

# Replace all occurrences
#lines <- gsub("[const std::vector<int>& state, const std::vectordouble>& p]", "(const std::vector<int>& state, const std::vectordouble>& p)", lines, fixed = TRUE)
#lines <- gsub("const std::vectordouble>& p", "const std::vector<double>& p", lines, fixed = TRUE)
#lines <- gsub("Eigen::VectorXd& state", "const Eigen::VectorXd& state", lines, fixed = TRUE)
lines <- gsub("[", "(", lines, fixed = TRUE)
lines <- gsub("]", ")", lines, fixed = TRUE)

# Write it back
writeLines(lines, "rates.txt")