library(Rcpp)

source("R/EvaluationFunctions.R")
source("R/FindClones.R")
sourceCpp("src/rcpp_hello_world.cpp")

data_path <- "/Users/seonghwanjun/data/simulation/binary4/"
for (case in 1:4) {
    sim_path <- paste(data_path, "/case", case, "/sim0", sep="")
    GetBSciteResults(sim_path, 10)
}
