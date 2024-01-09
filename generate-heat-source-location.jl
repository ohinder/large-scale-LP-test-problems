# https://en.wikipedia.org/wiki/Heat_equation
# **at equilbrium**
using Random
using ArgParse
import HDF5
using LinearAlgebra
using JuMP, SCS


function build_discretized_possion!(model, u, q, N::Int64)
    h = 1.0 / N
    @expression(model, grad_u_xx[i=1:N,j=1:N,k=1:N], (u[i+1,j,k] - 2 * u[i,j,k] + u[i-1,j,k]) / h^2);
    @expression(model, grad_u_yy[i=1:N,j=1:N,k=1:N], (u[i,j+1,k] - 2 * u[i,j,k] + u[i,j-1,k]) / h^2);
    @expression(model, grad_u_zz[i=1:N,j=1:N,k=1:N], (u[i,j,k+1] - 2 * u[i,j,k] + u[i,j,k-1]) / h^2);
    @constraint(model, grad_u_xx + grad_u_yy + grad_u_zz .== -q);

    # boundary conditions
    for j = 0:(N+1)
        for k = 0:(N+1)
            JuMP.fix(u[0,j,k], 0.0; force=true)
            JuMP.fix(u[N+1,j,k], 0.0; force=true)
        end
    end

    for i = 0:(N+1)
        for k = 0:(N+1)
            JuMP.fix(u[i,0,k], 0.0; force=true)
            JuMP.fix(u[i,N+1,k], 0.0; force=true)
        end
    end

    for i = 0:(N+1)
        for j = 0:(N+1)
            JuMP.fix(u[i,j,0], 0.0; force=true)
            JuMP.fix(u[i,j,N+1], 0.0; force=true)
        end
    end
end

function build_heat_source_detection_problem(
        num_source_locations::Int64,
        num_possible_source_locations::Int64,
        num_measurement_locations::Int64,
        grid_size::Int64
    )
    @assert num_source_locations < num_possible_source_locations
    @assert grid_size > 2
    @assert num_measurement_locations > num_source_locations
    # https://en.wikipedia.org/wiki/Inverse_problem

    heat_source_locations = rand(3, num_source_locations)
    heat_source_location_indexes = Int.(round.((grid_size - 1) * heat_source_locations)) .+ 1
    heat_flow_rate = rand(num_source_locations)

    q = zeros(grid_size, grid_size, grid_size)
    for location=axes(heat_source_location_indexes, 2)
        q[heat_source_location_indexes[:,location]...] = grid_size^3 * heat_flow_rate[location]
    end

    ##################
    # Compute u_true #
    ##################
    
    pde_model = Model(optimizer_with_attributes(SCS.Optimizer,
       "max_iters" => 100000,
       "eps" => 10^-11,
       "linear_solver" => SCS.IndirectSolver
    ))
    
    #pde_model = Model(HiGHS.Optimizer)
    @variable(pde_model, u[i=0:(grid_size+1),j=0:(grid_size+1),k=0:(grid_size+1)]);
    build_discretized_possion!(pde_model, u, q, grid_size)
    optimize!(pde_model)
    u_true = Array(value.(u))

    println("computed u_true")

    println("building inverse problem ...")
    measurement_locations = rand(3, num_measurement_locations)
    measurement_location_indexes = Int.(round.((grid_size - 1) * measurement_locations)) .+ 1
    
    candidate_locations = rand(3, num_possible_source_locations - num_source_locations)
    candidate_location_indexes = Int.(round.((grid_size - 1) * candidate_locations)) .+ 1
    

    ############################
    # Build optimization model #
    ############################
    model = Model(HiGHS.Optimizer)
    model = Model(optimizer_with_attributes(SCS.Optimizer,
       "max_iters" => 100000,
       "eps" => 10^-6,
       "linear_solver" => SCS.IndirectSolver
    ))
    @variable(model, u[i=0:(grid_size+1),j=0:(grid_size+1),k=0:(grid_size+1)]);
    @variable(model, 0.0 <= q[i=1:grid_size,j=1:grid_size,k=1:grid_size] <= 0.0);

    @objective(model, Min, sum(q));

    # second-order central difference
    build_discretized_possion!(model, u, q, grid_size)

    # q is unknown at possible heat source locations 
    for location=axes(heat_source_location_indexes, 2)
        location_indicies = heat_source_location_indexes[:,location] 
        JuMP.set_upper_bound(q[location_indicies...], grid_size^3)
    end
    for location=axes(candidate_location_indexes, 2)
        location_indicies = candidate_location_indexes[:,location]
        JuMP.set_upper_bound(q[location_indicies...], grid_size^3)
    end

    # u is known at the measurement locations 
    for location=axes(measurement_location_indexes,2)
        location_indicies = measurement_location_indexes[:,location] .+ 1
        JuMP.fix(u[location_indicies...], u_true[location_indicies...]; force=true)
    end
    println("built inverse problem")

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
        "--grid_size"
            help = "Number of grid planes per coordinate. The total number of grid points 
            will be the square of this number."
            arg_type = Int64
            default = 7
        "--ground_truth_file"
            help = "This is the file that will include the data for the ground truth value of u.
            This is a HDF5 file."
            arg_type = String
            required = true
        "--output_file"
            help = "This is the location that the mps file will be written to."
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    println("building model ...")

    Random.seed!(parsed_args["seed"])

    model, u_true = build_heat_source_detection_problem(
        parsed_args["num_source_locations"],
        parsed_args["num_possible_source_locations"],
        parsed_args["num_measurement_locations"],
        parsed_args["grid_size"]
    )

    HDF5.h5write(parsed_args["ground_truth_file"], "u_true", u_true)
    write_to_file(model, parsed_args["output_file"])
end

main()