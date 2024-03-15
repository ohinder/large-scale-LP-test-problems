# https://en.wikipedia.org/wiki/Heat_equation
# **at equilbrium**
using Random
using Dates
using ArgParse
import HDF5
using LinearAlgebra
using JuMP
using SparseArrays, IterativeSolvers
using Base
using HiGHS
include("utils.jl")

function solve_pde_linear_system(N::Int64, data::JuMP.LPMatrixData, pde_solve_tolerance::Float64)
    start_time = now()
    uvars = data.x_upper
    lvars = data.x_lower
    rhs = data.b_lower
    A = data.A
    varnames = [name(x) for x in data.variables]

    u_true_dict = Dict()
    keep_indicies = Vector{Int64}()
    for i in eachindex(lvars)
        lvar = lvars[i]
        if lvar == uvars[i]
            rhs -= A[:,i] * lvars[i]
            u_true_dict[varnames[i]] = lvar
        elseif lvar == -Inf && uvars[i] == Inf
            push!(keep_indicies, i)
        else
            error("Something was wrong with the model")
        end
    end
    rhs = Vector(rhs)
    A = A[:,keep_indicies]
    println("Generate ground truth system: ", now() - start_time)
    flush(stdout)

    start_time = now()
    u_tmp = minres(A, rhs, abstol=0.0, reltol = 1e-2 * pde_solve_tolerance)
    println("Solve ground truth system: ", now() - start_time)
    flush(stdout)

    if norm(A * u_tmp - rhs) > pde_solve_tolerance * norm(rhs)
        println("warning, didn't solve to desired accuracy:")
        @show norm(A * u_tmp - rhs) / norm(rhs)
    end

    for i in 1:length(keep_indicies)
        u_name = varnames[keep_indicies[i]]
        u_true_dict[u_name] = u_tmp[i]
    end


    u_true = zeros(N+2, N+2, N+2);
    for (key,value) in u_true_dict
       indicies = parse.(Int, split(key[3:end-1],","))
       u_true[indicies...] = value 
    end

    return u_true
end

function build_discretized_poisson!(model, u, q, N::Int64)
    h = 1.0 / N
    @expression(model, grad_u_xx[i=2:(N+1), j=2:(N+1), k=2:(N+1)], (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) / h^2)
    @expression(model, grad_u_yy[i=2:(N+1), j=2:(N+1), k=2:(N+1)], (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) / h^2)
    @expression(model, grad_u_zz[i=2:(N+1), j=2:(N+1), k=2:(N+1)], (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) / h^2)
    @constraint(model, grad_u_xx + grad_u_yy + grad_u_zz .== -q)

    # boundary conditions
    for j = 1:(N+2)
        for k = 1:(N+2)
            JuMP.fix(u[1, j, k], 0.0; force=true)
            JuMP.fix(u[N+2, j, k], 0.0; force=true)
        end
    end

    for i = 1:(N+2)
        for k = 1:(N+2)
            JuMP.fix(u[i, 1, k], 0.0; force=true)
            JuMP.fix(u[i, N+2, k], 0.0; force=true)
        end
    end

    for i = 1:(N+2)
        for j = 1:(N+2)
            JuMP.fix(u[i, j, 1], 0.0; force=true)
            JuMP.fix(u[i, j, N+2], 0.0; force=true)
        end
    end
end

function build_heat_source_detection_problem(
    num_source_locations::Int64,
    num_possible_source_locations::Int64,
    num_measurement_locations::Int64,
    grid_size::Int64,
    maximum_relative_measurement_error::Float64,
    optimize_model::Bool,
    pde_tolerance::Float64
)
    @assert num_source_locations < num_possible_source_locations
    @assert grid_size > 2
    @assert num_measurement_locations > num_source_locations
    # https://en.wikipedia.org/wiki/Inverse_problem

    start_time = now()
    heat_source_locations = rand(3, num_source_locations)
    heat_source_location_indexes = Int.(round.((grid_size - 1) * heat_source_locations)) .+ 2
    heat_flow_rate = rand(num_source_locations)

    q = zeros(grid_size, grid_size, grid_size)
    for location = axes(heat_source_location_indexes, 2)
        location_indicies = heat_source_location_indexes[:, location] .- 1
        q[location_indicies...] += grid_size^3 * heat_flow_rate[location]
    end
    println("Generate instance data: ", now() - start_time)
    flush(stdout)
    
    ##################
    # Compute u_true #
    ##################

    println("computing u_true")
    flush(stdout)

    start_time = now()
    pde_model = Model()
    @variable(pde_model, u[i=1:(grid_size+2), j=1:(grid_size+2), k=1:(grid_size+2)])
    build_discretized_poisson!(pde_model, u, q, grid_size)
    u_true = solve_pde_linear_system(grid_size, lp_matrix_data(pde_model), pde_tolerance)

    println("computed u_true: ", now() - start_time)
    flush(stdout)
    @show norm(u_true, 1)

    println("building inverse problem ...")
    flush(stdout)
    start_time = now()
    measurement_locations = rand(3, num_measurement_locations)
    measurement_location_indexes = Int.(round.((grid_size - 1) * measurement_locations)) .+ 2
    @assert minimum(measurement_location_indexes) >= 2
    @assert maximum(measurement_location_indexes) <= grid_size + 1

    candidate_locations = rand(3, num_possible_source_locations - num_source_locations)
    candidate_location_indexes = Int.(round.((grid_size - 1) * candidate_locations)) .+ 2
    @assert maximum(candidate_location_indexes) <= grid_size + 1
    @assert minimum(candidate_location_indexes) >= 2
    println("Built inverse problem data: ", now() - start_time)
    flush(stdout)


    ############################
    # Build optimization model #
    ############################
    println("building model ...")
    flush(stdout)

    begin
        model_start_time = now()
        if optimize_model
            model = Model(HiGHS.Optimizer)
        else 
            model = Model()
        end
	start_time = now()
        @variable(model, u[i=1:(grid_size+2), j=1:(grid_size+2), k=1:(grid_size+2)])
        @variable(model, 0.0 <= q[i=1:grid_size, j=1:grid_size, k=1:grid_size] <= 0.0)
	println("Create variables: ", now() - start_time)
	flush(stdout)

        start_time = now()
        @objective(model, Min, sum(q))
	println("Set objective: ", now() - start_time)
	flush(stdout)

	start_time = now()
        # second-order central difference
        build_discretized_poisson!(model, u, q, grid_size)
	println("Create discretized poisson constraints: ", now() - start_time)
	flush(stdout)

        # q is unknown at possible heat source locations
        for location = axes(heat_source_location_indexes, 2)
            location_indicies = heat_source_location_indexes[:, location] .- 1
            if has_upper_bound(q[location_indicies...])
                JuMP.delete_upper_bound(q[location_indicies...])
            end
        end
        for location = axes(candidate_location_indexes, 2)
            location_indicies = candidate_location_indexes[:, location] .- 1
            if has_upper_bound(q[location_indicies...])
                JuMP.delete_upper_bound(q[location_indicies...])
            end
        end

        # u is known at the measurement locations
        for location = axes(measurement_location_indexes, 2)
            location_indicies = measurement_location_indexes[:, location]
            # Allow a small amount of errors in measurments
            # to ensure that errors in the calculation of u_true do not
            # make the model infeasible.
            maximum_u_value = u_true[location_indicies...] * (1 + maximum_relative_measurement_error)
            minimum_u_value = u_true[location_indicies...] / (1 + maximum_relative_measurement_error)

            JuMP.set_lower_bound(u[location_indicies...], minimum_u_value)
            JuMP.set_upper_bound(u[location_indicies...], maximum_u_value)
        end
        println("built inverse problem: ", now() - model_start_time)
	flush(stdout)
    end

    if optimize_model
        println("optimizing model ...")
        optimize!(model)
        u_opt = value.(u)
        println("")
        println("How well did we recover the ground truth:")
        @show norm(u_opt - u_true) / norm(u_true)
    end
    return model, u_true
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        help = "The random seed used to generate the instance."
        arg_type = Int64
        default = 1
        "--num_source_locations"
        arg_type = Int64
        default = 2
        "--num_possible_source_locations"
        arg_type = Int64
        default = 50
        "--num_measurement_locations"
        arg_type = Int64
        default = 20
        "--maximum_relative_measurement_error"
        arg_type = Float64
        default = 1e-6
        "--grid_size"
        help = "Number of grid planes per coordinate. The total number of grid points
        will be the square of this number."
        arg_type = Int64
        default = 7
        "--pde_solve_tolerance"
        arg_type = Float64
        default = 1e-8
        "--optimize_model"
        arg_type = Bool
        default = false
        "--ground_truth_file"
        help = "This is the file that will include the data for the ground truth value of u.
        This is a HDF5 file."
        arg_type = String
        required = true
        "--output_file"
        help = "This is the location that the mps file will be written to."
        required = true
        "--rescale_model"
        help = RESCALE_MODEL_HELP_DESCRIPTION # this is defined in utils.jl
        arg_type = Bool
        default = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg, val) in parsed_args
        println("  $arg  =>  $val")
    end

    Random.seed!(parsed_args["seed"])

    @time "Build model" model, u_true = build_heat_source_detection_problem(
        parsed_args["num_source_locations"],
        parsed_args["num_possible_source_locations"],
        parsed_args["num_measurement_locations"],
        parsed_args["grid_size"],
        parsed_args["maximum_relative_measurement_error"],
        parsed_args["optimize_model"],
        parsed_args["pde_solve_tolerance"]
    )
    flush(stdout)

    if isfile(parsed_args["ground_truth_file"])
        rm(parsed_args["ground_truth_file"])
    end
    @time "Write ground truth" HDF5.h5write(parsed_args["ground_truth_file"], "u_true", u_true)
    flush(stdout)

    if parsed_args["rescale_model"]
        println("rescaling model ...")
	flush(stdout)
        @time "Rescale model" model = rescale_instance(lp_matrix_data(model))
	flush(stdout)
    end

    @time "Write model" write_to_file(model, parsed_args["output_file"])
    flush(stdout)
end

main()
