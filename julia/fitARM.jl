include("armoptimize.jl")

df = CSV.read("raw_data/anes_timeseries_cdf_csv_20220916/anes_timeseries_cdf_csv_20220916.csv", DataFrame);
filtered_df = df[(df[:, "VCF0302"] .== "1") .| (df[:, "VCF0302"] .== "5"), :];

var_map = SortedDict{String, String}(
    "VCF0834" => "Women should be in the home?", # 7
    "VCF0838" => "Should abortion be legal?", #16
    "VCF0876a" => "Laws shouldn't protect homosexuality?", #9
    "VCF0879" => "Decrease number of immigrants?", #9
    "VCF9043" => "Should school prayer be allowed?", #7
    "VCF9237" => "Oppose death penalty?", #5
    "VCF9238" => "Should be easier to buy guns?", #6
    "VCF0830" => "Blacks shouldn't have aid?", #17
    "VCF0867a" => "Oppose Affirmative Action?", #12
    "VCF0815" => "Segregation", #6
    "VCF0817" => "There shouldn't be school busing for integration?", #5
    "VCF0842" => "Environmental Regulation is too burdensome?" #4
);

full_results = Dict{String, Any}()
for variable in keys(var_map)

    res = get_series(filtered_df, variable)
    full_results[variable] = res

end

println("number threads: ", Threads.nthreads())

const E = [i for i in .05:.05:.95]
const T = [i for i in .05:.05:.95]
const R = [i for i in .05:.05:.95]

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, false, false)

js = JSON.json(optimals)

open("outputs/base_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, true, false)

js = JSON.json(optimals)

open("outputs/id_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, false, true)

js = JSON.json(optimals)

open("outputs/h_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, true, true)

js = JSON.json(optimals)

open("outputs/id_h_optimals.json", "w") do f
    write(f, js)
end

# julia --threads 10 fitARM.jl