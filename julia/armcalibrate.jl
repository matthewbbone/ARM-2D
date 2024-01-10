
function simulate(N::Int64, n_agents::Int64, cycles::Int64, affect_ad::Bool, affect_h::Bool, locked::Vector{Float64}=[-1., -1., -1.])
    
    results::Vector{Any} = zeros(N)
    lk = ReentrantLock()

    Threads.@threads for thread=1:N
        E = locked[1] == -1. ? rand(Beta(1,1)) : locked[1] 
        T = locked[2] == -1. ? rand(Beta(1,1)) : locked[2] 
        R = locked[3] == -1. ? rand(Beta(1,1)) : locked[3] 
        agents = [Agent(i, limit(randn() / 5 + .5, "norm"), E, T, R, rand(['A','B'])) for i=1:n_agents]
        res = [runARM(cycles, agents, affect_ad, affect_h)]
        lock(lk) do
            results[thread] = [E, T, R, res]
        end
    end

    return results

end 

function distance(sample, truth)
    return sum((sample[1][4] - truth[1][4]).^2)
end

function plot_posterior(posterior::Vector{Float64}, truth::Float64)
    p = plot(histogram(x=posterior, nbinsx=100, range_x=[0,1]), Layout(xaxis_range=[0, 1], width=500))
    
    add_shape!(p, line(
        x0=mean(posterior), y0=0,
        x1=mean(posterior), y1=length(posterior),
        line=attr(color="Red", width=3),
    ))

    add_shape!(p, line(
        x0=truth, y0=0,
        x1=truth, y1=length(posterior),
        line=attr(color="Green", width=3),
    ))

    return p
end

function abc(truth::Vector{Any}, n_samples::Int64, threshold::Float64, affect_ad::Bool, affect_h::Bool, locked::Vector{Float64})
    
    posterior::Vector{Float64} = []
    for i in 1:n_samples
        dist = threshold + 1
        theta = 0
        while dist > threshold
            theta = rand(Beta(1,1))
            sample = getAvgResults(simulate(1, 100, 2500, affect_ad, affect_h, [theta, .25, .25]))
            dist = distance(truth, sample)
        end
        push!(posterior, theta)
    end
    return posterior
end

# https://link.springer.com/article/10.1007/s11222-011-9296-2
function metropolis_hastings(truth::Vector{Any}, n_samples::Int64, n_iters::Int64, affect_ad::Bool, affect_h::Bool, locked::Vector{Float64})

    thetas::Vector{Float64} = []
    distances::Vector{Float64} = []
    ctr = 0

    for i in 1:n_samples

        theta = rand(Beta(1,1))
        sample = getAvgResults(simulate(1, 100, 2500, affect_ad, affect_h, [theta, .25, .25]))
        dist = distance(truth, sample)
        
        push!(thetas, theta)
        push!(distances, dist)

    end

    theta_sets = [thetas]

    for i in 1:n_iters

        old_kde = kde(thetas, boundary=(0,1), bandwidth=.1)

        z = sample(old_kde.x, Weights(old_kde.density))
        z_sample = getAvgResults(simulate(1, 100, 2500, affect_ad, affect_h, [z, .25, .25]))
        z_dist = distance(truth, z_sample)

        n_index = rand(1:length(thetas))
        n = thetas[n_index]

        prop_thetas = copy(thetas)
        prop_thetas[n_index] = z
        prop_kde = kde(prop_thetas, boundary=(0,1), bandwidth=.1)

        z_kde = pdf(old_kde, z)
        n_kde = pdf(prop_kde, n)

        proxy_dist = (1 ./ distances) / sum(1 ./ distances)

        prop_distances = copy(distances)
        prop_distances[n_index] = z_dist
        prop_proxy_dist = (1 ./ prop_distances) / sum(1 ./ prop_distances)

        switch_prob = (prop_proxy_dist[n_index] * n_kde) / (proxy_dist[n_index] * z_kde)
        # switch_prob = prop_proxy_dist[n_index] / proxy_dist[n_index]

        if rand(Beta(1,1)) < min(1, switch_prob)
            ctr = ctr + 1
            thetas = prop_thetas
            distances = prop_distances
        end

        push!(theta_sets, thetas)

    end

    posterior = theta_sets[1]
    for i in 2:length(theta_sets)
        posterior = posterior .+ theta_sets[i]
    end

    return thetas
end