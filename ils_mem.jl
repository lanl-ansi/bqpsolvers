#!/usr/bin/env julia

### Note ###
#
# code from the blog post at, https://invenia.github.io/blog/2019/09/11/memristors/
#
# tested on Julia v1.2
# test with: include("ils_mem.jl"); main(Dict("random-seed" => nothing, "max-steps"=>100, "show-solution" => true, "runtime-limit"=>1.0, "input-file" => "---"))
#

using ArgParse
using JSON

using Random
using LinearAlgebra

function main(parsed_args)
    if !isnothing(parsed_args["random-seed"])
        Random.seed!(parsed_args["random-seed"])
    end

    #read the problem
    file = open(parsed_args["input-file"])
    data = JSON.parse(file)

    if data["variable_domain"] != "boolean"
        error("only boolean domains are supported. Given $(data["variable_domain"])")
    end

    idx_to_var = data["variable_ids"]
    var_to_idx = Dict(var => idx for (idx,var) in enumerate(data["variable_ids"]))
    n = length(idx_to_var)

    #println(idx_to_var)
    println("problem size: $(n), $(length(data["quadratic_terms"]))")

    Q = zeros(n, n);
    for qt in data["quadratic_terms"]
        i = var_to_idx[qt["id_head"]]
        j = var_to_idx[qt["id_tail"]]
        c = qt["coeff"]
        Q[i, j] = c/2
        Q[j, i] = c/2
    end

    h = [0.0 for i in 1:n]
    for lt in data["linear_terms"]
        i = var_to_idx[lt["id"]]
        c = lt["coeff"]

        h[i] = c
    end

    #optional regularisation term
    for i in 1:n
        Q[i,i] = 1e-4
    end

    # complication run
    memristive_opt(h, Q, 0.1, [0.5 for i in 1:n], total_time=1)


    #v = [0 for i in 1:n]
    #check_encoding(h, Q, v)

    #v = [1 for i in 1:n]
    #check_encoding(h, Q, v, data)

    #v = [mod(i,2) for i in 1:n]
    #check_encoding(h, Q, v, data)

    # a tunable parameter
    p = 0.1

    restarts = 0
    best_assignment = Dict(idx_to_var[i] => 0.0 for i in 1:n)
    best_energy = calc_energy(data, best_assignment)
    print("($(best_energy))")

    time_start = time()
    while time() - time_start < parsed_args["runtime-limit"]
        v = [Random.rand() for i in 1:n]

        weights_final, energies = memristive_opt(h, Q, p, v, total_time=parsed_args["max-steps"])

        #println("energies trace: $energies")
        #println("weights: $weights_final")

        assignment = Dict(idx_to_var[i] => weights_final[i] for i in 1:n)
        energy = calc_energy(data, assignment)

        if energy < best_energy
            best_energy = energy
            best_assignment = assignment
            print("($(best_energy))")
        end

        restarts += 1
        print("R")
    end
    println()

    time_elapsed = time() - time_start

    println("restarts: $restarts")

    println("final energy eval: $best_energy")

    if parsed_args["show-solution"]
        println()
        for (i,v) in sort(collect(best_assignment), by=first)
            println("$(i) - $(v)")
        end
    end


    nodes = length(data["variable_ids"])
    edges = length(data["quadratic_terms"])

    scale = data["scale"]
    offset = data["offset"]
    lt_lb = -sum(abs(lt["coeff"]) for lt in data["linear_terms"])/scale
    qt_lb = -sum(abs(qt["coeff"]) for qt in data["quadratic_terms"])/scale
    lower_bound = lt_lb+qt_lb

    best_objective = best_energy
    best_nodes = restarts
    best_runtime = time_elapsed
    scaled_objective = scale*(best_objective+offset)
    scaled_lower_bound = scale*(lower_bound+offset)

    println("BQP_DATA, $(nodes), $(edges), $(scaled_objective), $(scaled_lower_bound), $(best_objective), $(lower_bound), $(best_runtime), $(0), $(best_nodes)")
end


"""
    memristive_opt(
        expected_returns::Vector{Float64},
        Σ::Matrix{Float64},
        p::Float64,
        weights_init::Vector{<:Real};
        α=0.1,
        β=10,
        δt=0.1,
        total_time=3000,
    )

Execute optimisation via the heuristic "memristive" equation in order to find the
optimal portfolio composition, considering an asset covariance matrix `Σ` and a
risk parameter `p`. `reg` represents the regularisation constant for `Σ`, `α` and
`β` parametrise the memristor state (see Equation (2)),
`δt` is the size of the time step for the dynamical updates and `total_time` is
the number of time steps for which the dynamics will be run.
"""
function memristive_opt(
    expected_returns::Vector{Float64},
    Σ::Matrix{Float64},
    p::Float64,
    weights_init::Vector{<:Real};
    α=0.1,
    β=10,
    δt=0.1,
    total_time=3000,
)
    n = size(weights_init, 1)

    weights_series = Matrix{Float64}(undef, n, total_time)
    weights_series[:, 1] = weights_init
    energies = Vector{Float64}(undef, total_time)
    energies[1] = energy(expected_returns, Σ, p, weights_series[:, 1])

    # Compute resistance change ratio
    ξ = p / 2α

    #println(p, " ", α, " ", β, " ", ξ)
    #display(Σ)
    #S = β * inv(Σ) * (α/2 * ones(n, 1) + (p/2 + α * ξ/3) * diag(Σ) - expected_returns)
    #println(S)

    # Compute Σ times applied voltages matrix
    ΣS = β * (α/2 * ones(n, 1) + (p/2 + α * ξ/3) * diag(Σ) - expected_returns)

    #println(inv(Σ)*ΣS)
    prev_update = zeros(n)
    t = 1
    while t < total_time
    #for t in 1:total_time-1
        update = δt * (α * weights_series[:, t] - 1/β * (I(n) + ξ *
            Σ * Diagonal(weights_series[:, t])) \ ΣS)

        #println(norm(update .- prev_update))
        if norm(update .- prev_update) <= 1e-6
            break # solver converged
        else
            #println(norm(update .- prev_update), " ", update)
        end
        prev_update = update

        weights_series[:, t+1] = weights_series[:, t] + update

        weights_series[weights_series[:, t+1] .> 1, t+1] .= 1.0
        weights_series[weights_series[:, t+1] .< 1, t+1] .= 0.0

        energies[t + 1] = energy(expected_returns, Σ, p, weights_series[:, t+1])
        #println(weights_series[:, t+1])
        #println(energies[t + 1])
        t += 1
    end

    #println(t)
    weights_final = weights_series[:, t]

    #println(weights_final)
    #println(energies)

    return weights_final, energies
end


"""
    energy(
        expected_returns::Vector{Float64},
        Σ::Matrix{Float64},
        p::Float64,
        weights::Vector{Float64},
    )

Return minus the expected return corrected by the variance of the portfolio,
according to the Markowitz approach. `Σ` represents the covariance of the assets,
`p` controls the risk tolerance and `weights` represent the (here binary)
portfolio composition.
"""
function energy(
    expected_returns::Vector{Float64},
    Σ::Matrix{Float64},
    p::Float64,
    weights::Vector{<:Real},
)
    -dot(expected_returns, weights) + p/2 * weights' * Σ * weights
end


"checks that qubo evalution matches energy function evaluation"
function check_encoding(h, Q, assignment, data)
    eval_q = assignment' * Q * assignment - dot(h, assignment)
    eval_energy = energy(h, Q, 2.0, assignment)
    #println(eval_q, " ", eval_energy)

    #eval_data = calc_energy(data, Dict(data["variable_ids"][i] => v for (i,v) in enumerate(assignment)))
    #println(eval_q, " ", eval_energy, " ", eval_data)
    @assert(isapprox(eval_q, eval_energy))
    #@assert(isapprox(eval_energy, eval_data))
end


"evaluate the enery function given a variable assignment"
function calc_energy(data, assignment)::Float64
    energy = 0.0
    for qt in data["quadratic_terms"]
        i = qt["id_head"]
        j = qt["id_tail"]
        c = qt["coeff"]
        energy += c * assignment[i] * assignment[j]
    end

    for lt in data["linear_terms"]
        i = lt["id"]
        c = lt["coeff"]
        energy += c * assignment[i]
    end

    #return data["scale"]*(energy+data["offset"])
    return energy
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input-file", "-f"
            help = "the data file to operate on (.json)"
            required = true
        "--runtime-limit", "-t"
            help = "puts a time limit on the solver"
            arg_type = Float64
            default = 10.0
        "--show-solution", "-s"
            help = "print the solution"
            action = :store_true
        "--random-seed"
            help = "fixes the random number generator seed"
            arg_type = Int
        "--max-steps"
            help = "set the maximum number of steps in the memristor evolution"
            arg_type = Int
            default = 100
    end

    return parse_args(s)
end

#    parser.add_argument('-ss', '--show-solution', help='print the solution', action='store_true', default=False)
#    parser.add_argument('-rtl', '--runtime-limit', help='gurobi runtime limit (sec.)', type=float)


if isinteractive() == false
    main(parse_commandline())
end
