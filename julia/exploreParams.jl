# This file can be run with the command "julia --threads 10 armtest.jl"
# It runs a parameter grid search of the Attraction-Repulsion model presented
# in Axelrod et al. (https://www.pnas.org/doi/abs/10.1073/pnas.2102139118).
# A few modifications were made for this model to reflect affective polarization
# instead of ideological polarization. The development and exploration of this
# code can be reviewed here: https://github.com/matthewbbone/ARM-2D/tree/main/julia
# I've separately added the working paper that is reporting on the dynamics of this
# model relative Axelrod et al.'s original model.

include("ARM.jl")

println("number threads: ", Threads.nthreads())

const E = [i for i in .05:.05:.95]
const T = [i for i in .05:.05:.95]
const R = [i for i in .05:.05:.95]

@time results = constantSuite(30, 500, 100, E, T, R, false, true, true);

# avg_results = getAvgResults(results);

js = JSON.json(results)

open("h_affect_results.json", "w") do f
    write(f, js)
end

# julia --threads 10 exploreParams.jl 