using PlotlyJS, Statistics, BenchmarkTools, DataFrames, JSON, StatsBase, Random, Distributions, KernelDensity, ProgressBars, CSV, DataFramesMeta, PlotlyBase, DataStructures
include("armtest.jl")

function get_series(df, variable)

    df = deepcopy(df)

    filter!(row -> row[variable] != "9", df)
    filter!(row -> row[variable] != "0", df)
    filter!(row -> row[variable] != " ", df)
    filter!(row -> row[variable] != "-8", df)
    filter!(row -> row[variable] != "-9", df)
    filter!(row -> row[variable] != "8", df)
    filter!(row -> row[variable] != "9", df)

    if (variable .== "VCF0876a") .| (variable .== "VCF0867a")
        filter!(row -> row[variable] != "7", df)
    end

    transform!(df, variable => (x -> parse.(Int, x)) => variable)
    transform!(df, variable => (x -> (x .- minimum(x)) ./ (maximum(x) - minimum(x))) => variable)

    df = select(df, "VCF0004", "VCF0302", variable)

    years = unique(df[!,"VCF0004"])

    stats = Dict{String,Any}(
        "mean" => [],
        "variance" => [],
        "range" => [],
        "diffs" => []
    )

    dem_stats = deepcopy(stats)
    rep_stats = deepcopy(stats)

    stats["perc_dem"] = size(filter(row -> row["VCF0302"] == "1", df), 1) / size(df, 1)
    stats["perc_rep"] = size(filter(row -> row["VCF0302"] == "5", df), 1) / size(df, 1)

    @assert stats["perc_dem"] + stats["perc_rep"] == 1

    for t in eachindex(years)

        temp_df = filter(row -> row["VCF0004"] == years[t], df)
        rep_df = filter(row -> row["VCF0302"] == "1", temp_df)
        dem_df = filter(row -> row["VCF0302"] == "5", temp_df)

        push!(stats["mean"], mean(temp_df[!, variable]))
        push!(stats["variance"], var(temp_df[!, variable]))
        push!(stats["range"], temp_df[!, variable])

        push!(dem_stats["mean"], mean(dem_df[!, variable]))
        push!(dem_stats["variance"], var(dem_df[!, variable]))
        push!(dem_stats["range"], dem_df[!, variable])

        push!(rep_stats["mean"], mean(rep_df[!, variable]))
        push!(rep_stats["variance"], var(rep_df[!, variable]))
        push!(rep_stats["range"], rep_df[!, variable])

        if t == 1
            push!(stats["diffs"], 0)
            push!(dem_stats["diffs"], 0)
            push!(rep_stats["diffs"], 0)
        else
            push!(stats["diffs"], stats["mean"][t] - stats["mean"][t-1])
            push!(dem_stats["diffs"], dem_stats["mean"][t] - dem_stats["mean"][t-1])
            push!(rep_stats["diffs"], rep_stats["mean"][t] - rep_stats["mean"][t-1])
        end
        
    end

    results = Dict(
        "variable" => variable,
        "years" => years,
        "stats" => stats,
        "dem_stats" => dem_stats,
        "rep_stats" => rep_stats
    )

    return results

end 

function display_series(results, var_map)

    p = plot(Layout(title=var_map[results["variable"]], yaxis=attr(range=[0,1])))

    mean_trace = scatter(
        x=results["years"], 
        y=results["stats"]["mean"], mode="lines", name="mean",
        marker_color = "black"
    )

    dem_trace = scatter(
        x=results["years"], 
        y=results["dem_stats"]["mean"], mode="lines", name="dem mean",
        marker_color = "blue"
    )

    rep_trace = scatter(
        x=results["years"], 
        y=results["rep_stats"]["mean"], mode="lines", name="rep mean",
        marker_color = "red"
    )

    return [mean_trace, dem_trace, rep_trace]

end

function compare_results(
    E::Float64, 
    T::Float64, 
    R::Float64, 
    affect_ad::Bool, 
    affect_h::Bool, 
    empirics::Dict{String, Any},
    var_map::SortedDict{String, String},
    pop_size::Int64)

    losses = Dict{String, Float64}("total" => 0.)
    for variable in keys(var_map)

        n_dems = floor(empirics[variable]["stats"]["perc_dem"] * pop_size)
        mean_dems = empirics[variable]["dem_stats"]["mean"][1]
        var_dems = empirics[variable]["dem_stats"]["variance"][1]
        mean_reps = empirics[variable]["rep_stats"]["mean"][1]
        var_reps = empirics[variable]["rep_stats"]["variance"][1]

        # one cycle is a month
        n_steps = (empirics[variable]["years"][end] - empirics[variable]["years"][1]) * 12 + 1
        
        agents::Vector{Agent} = []
        for i in 1:pop_size

            if i < n_dems
                push!(agents, Agent(i, limit(randn() * var_reps^.5  + mean_reps, "norm"), E, T, R, 'A'))
            else
                push!(agents, Agent(i, limit(randn() * var_dems^.5 + mean_dems, "norm"), E, T, R, 'B'))
            end

        end

        simulation = runARM(n_steps, agents, affect_ad, affect_h)
        res_mat = hcat(simulation...)

        rep_group = []
        dem_group = []
        total_group = []

        for i in 1:size(res_mat,1)
            if agents[i].group == 'A'
                push!(rep_group, res_mat[i,:])
            else
                push!(dem_group, res_mat[i,:])
            end
            push!(total_group, res_mat[i,:])
        end

        rep_mean = mean(rep_group, dims=1)[1]
        rep_var = var.(eachcol(rep_group))[1]

        dem_mean = mean(dem_group, dims=1)[1]
        dem_var =var.(eachcol(dem_group))[1]

        total_mean = mean(total_group, dims=1)[1]
        total_var = var.(eachcol(total_group))[1]

        losses[variable] = 0.
        old_idx = 1
        for t in eachindex(empirics[variable]["years"][2:end])

            idx = (empirics[variable]["years"][t] - empirics[variable]["years"][1]) * 12 + 1
            
            # losses[variable] += abs(rep_mean[idx] - empirics[variable]["rep_stats"]["mean"][t])
            # losses[variable] += abs(rep_var[idx] - empirics[variable]["rep_stats"]["variance"][t])
            # losses[variable] += abs((rep_mean[idx] - rep_mean[old_idx]) - empirics[variable]["rep_stats"]["diffs"][t])

            # losses[variable] += abs(dem_mean[idx] - empirics[variable]["dem_stats"]["mean"][t])
            # losses[variable] += abs(dem_var[idx] - empirics[variable]["dem_stats"]["variance"][t])
            # losses[variable] += abs((dem_mean[idx] - dem_mean[old_idx]) - empirics[variable]["dem_stats"]["diffs"][t])

            # losses[variable] += abs(total_mean[idx] - empirics[variable]["stats"]["mean"][t])
            # losses[variable] += abs(total_var[idx] - empirics[variable]["stats"]["variance"][t])
            # losses[variable] += abs((total_mean[idx] - total_mean[old_idx]) - empirics[variable]["stats"]["diffs"][t])

            old_idx = idx

            rep_options = unique(empirics[variable]["rep_stats"]["range"][1])
            rep_sims = countmap([rep_options[argmin(abs.(rv[idx] .- rep_options))] for rv in rep_group])
            rep_sims = Dict(k => float(v) for (k, v) in rep_sims)
            rep_emp = countmap(empirics[variable]["rep_stats"]["range"][1])
            rep_emp = Dict(k => float(v) for (k, v) in rep_emp)

            dem_options = unique(empirics[variable]["dem_stats"]["range"][1])
            dem_sims = countmap([dem_options[argmin(abs.(rv[idx] .- dem_options))] for rv in dem_group])
            dem_sims = Dict(k => float(v) for (k, v) in dem_sims)
            dem_emp = countmap(empirics[variable]["dem_stats"]["range"][1])
            dem_emp = Dict(k => float(v) for (k, v) in dem_emp)

            for k in keys(rep_emp)
                if k in keys(rep_sims)
                    sim = rep_sims[k] / (pop_size - n_dems)
                    emp = rep_emp[k] / length(empirics[variable]["rep_stats"]["range"][1])
                    losses[variable] += abs(emp - sim)
                end
                if k in keys(dem_sims)
                    sim = dem_sims[k] / n_dems
                    emp = dem_emp[k] / length(empirics[variable]["dem_stats"]["range"][1])
                    losses[variable] += abs(emp - sim)
                end
            end
            
        end
        losses[variable] = (losses[variable] / length(empirics[variable]["years"][2:end]))
        losses["total"] +=  losses[variable]

    end

    losses["total"] /=  length(var_map)

    return losses

end


function get_optimals(
    trials::Int64,
    n_agents::Int64, 
    empirics::Dict{String, Any},
    var_map::SortedDict{String, String},
    E::Vector{Float64}, 
    T::Vector{Float64}, 
    R::Vector{Float64}, 
    affect_ad::Bool, 
    affect_h::Bool)

    results::Vector{Any} = zeros(length(E) * length(T) * length(R))
    lk = ReentrantLock()
    
    for e in ProgressBar(eachindex(E))
        Threads.@threads for t in eachindex(T)
            for r in eachindex(R)
                res = compare_results(E[e], T[t], R[r], affect_ad, affect_h, empirics, var_map, n_agents)
                for tr=2:trials
                    temp_res = compare_results(E[e], T[t], R[r], affect_ad, affect_h, empirics, var_map, n_agents)
                    for k in keys(res)
                        res[k] += temp_res[k]
                    end
                end

                for k in keys(res)
                    res[k] /= trials
                end

                lock(lk) do
                    index = (e - 1) * length(T) * length(R) + (t - 1) * length(R) + r
                    results[index] = [E[e], T[t], R[r], res]
                end
            end
        end
    end

    return results

end

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

plots = []
p = make_subplots(rows=floor(Int, length(var_map) / 4), cols=4)

global ctr = 0
full_results = Dict{String, Any}()
for variable in keys(var_map)

    res = get_series(filtered_df, variable)
    full_results[variable] = res
    traces = display_series(res, var_map)
    
    for t in traces
        add_trace!(p, t, row=floor(Int, ctr / 4) + 1, col=ctr % 4 + 1)
    end

    global ctr += 1

end

println("number threads: ", Threads.nthreads())

const E = [i for i in .05:.05:.95]
const T = [i for i in .05:.05:.95]
const R = [i for i in .05:.05:.95]

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, false, false)

js = JSON.json(optimals)

open("base_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, true, false)

js = JSON.json(optimals)

open("id_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, false, true)

js = JSON.json(optimals)

open("h_optimals.json", "w") do f
    write(f, js)
end

@time optimals = get_optimals(10, 100, full_results, var_map, E, T, R, true, true)

js = JSON.json(optimals)

open("id_h_optimals.json", "w") do f
    write(f, js)
end

# julia --threads 10 armoptimize.jl