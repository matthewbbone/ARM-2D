using PlotlyJS
using Statistics
using BenchmarkTools 
using DataFrames
using JSON
using StatsBase
using Random
using Distributions 
using KernelDensity
using ProgressBars
using CSV
using DataFramesMeta
using PlotlyBase 
using DataStructures

include("ARM.jl")

function ideological_reaction(T)

    d_values = -1:0.01:1  # Range of d values

    p = make_subplots(
      rows=3, 
      cols=1, 
      shared_xaxes=true,
      shared_yaxes=true,
      row_titles=["original", "ingroup", "outgroup"]
    )
    
    # Create a plot
    for r in [.05, .15, .25, .35, .45, .55, .65, .75, .85, .95]
        
        i_values = [abs(d) <= T ? r*d : -r*d for d in d_values]
        
        trace = scatter(x=d_values, y=i_values, mode="lines")
        add_trace!(p, trace, row=1, col=1)
        
    end

    for r in [.05, .15, .25, .35, .45, .55, .65, .75, .85, .95]
        
        i_in = [r * abs(d) / d for d in d_values]

        trace = scatter(x=d_values, y=i_in, mode="lines")
        add_trace!(p, trace, row=2, col=1)
        
    end

    for r in [.05, .15, .25, .35, .45, .55, .65, .75, .85, .95]
        
        i_out = [- r * (abs(d) * (1-T)) / d for d in d_values]

        trace = scatter(x=d_values, y=i_out, mode="lines")
        add_trace!(p, trace, row=3, col=1)
        
    end

    relayout!(p,
        title="Reaction Distribution (T=$(T))", 
        showlegend=false
    )
    
    return p
end

function ideological_exposure()

    d_values = 0:0.01:1  # Range of d values
    p = plot()
    
    # Create a plot
    for e in [.05, .15, .25, .35, .45, .55, .65, .75, .85, .95]
        p_values = [(.5)^(d/e) for d in d_values]
        trace = scatter(;x=d_values, y=p_values, mode="lines", name="E=$(e)")
        addtraces!(p, trace)
    end

    line_trace = scatter(;x=[1/3, 1/3], y=[0, 1], mode="lines", line=attr(color="black"), name="outgroup exposure")
    addtraces!(p, line_trace)

    relayout!(p,
        title="Exposure Distribution", 
        xaxis_title="d", 
        yaxis_title="p"
    )
    return p
end

function cohens_d(group1, group2)
    n1, n2 = length(group1), length(group2)
    var1, var2 = var(group1), var(group2)
    pooled_var = ((n1 - 1)*var1 + (n2 - 1)*var2) / (n1 + n2 - 2)
    d = abs(mean(group1) - mean(group2)) / sqrt(pooled_var)
    return d
end

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
    stats["cohens"] = []

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

        push!(stats["cohens"], cohens_d(rep_df[!, variable], dem_df[!, variable]))
        
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

function optimalMap(optimals::Vector{Any}, variable::String, E::Float64, T::Float64, R::Float64, range1::Vector{Float64}, range2::Vector{Float64}, inorder::Bool)

    heat_dat = Array{Float64, 2}(undef, 0, length(range2))
    for i in range2
        row::Vector{Float64} = []
        for j in range1
            if inorder
                if E >= 0 p = searchRes(optimals, E, i, j)[4][variable]
                elseif T >= 0 p = searchRes(optimals, i, T, j)[4][variable]
                else p = searchRes(optimals, i, j, R)[4][variable]
                end
            else
                if E >= 0 p = searchRes(optimals, E, j, i)[4][variable]
                elseif T >= 0 p = searchRes(optimals, j, T, i)[4][variable]
                else p = searchRes(optimals, j, i, R)[4][variable]
                end
            end
            push!(row, p)
        end
        heat_dat = vcat(heat_dat, reshape(row, 1, length(row)))
    end

    return heatmap(z=heat_dat, x=range1, y=range2, zmin=1.3, zmax=1.7, colorscale=[[0, "rgb(0,128,0)"], [1, "rgb(255,0,0)"]])
end

