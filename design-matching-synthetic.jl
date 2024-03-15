using JuMP
import HiGHS
using Random
using Dates
using LinearAlgebra
using SparseArrays
using ArgParse
using NearestNeighbors
include("utils.jl")

# aim for problems with > 100 million nonzeros
# create a suite a smaller test problems with ~1 million nonzeros

# TODO(ohinder): consider using https://archive.ics.uci.edu/dataset/31/covertype to create 
# a real instance.

function build_synthetic_design_matching_problem(
    epsilon::Float64,
    num_treatment_samples::Int64,
    num_control_samples::Int64,
    num_covariates::Int64,
    num_edges_per_treatment::Int64,
    control_shift_magnitude::Float64,
    optimize_model::Bool
)
    start_time = now()
    n = num_treatment_samples
    m = num_control_samples
    @assert m > n
    @assert num_edges_per_treatment <= m
    d = num_covariates
    A = Random.randn(n, d); 
    shift = control_shift_magnitude * randn(d)
    B = Random.randn(m, d) + repeat(shift', m); # add shift

    

    # select only the closest `num_edges_per_treatment' edges per treatment
    kdtree = KDTree(B')
    edges = Vector{Tuple{Int64, Int64}}(undef, n * num_edges_per_treatment);
    col = 1
    for i = 1:n
        idxs, _ = knn(kdtree, A[i,:], num_edges_per_treatment)
        for j in idxs
            edges[col] = (i,j)
            col += 1
        end
    end

    c = Vector{Float64}(undef, length(edges));
    for col=1:length(edges) #(i,j) in edges
        (i,j) = edges[col]
        c[col] = norm(A[i,:] - B[j,:])
    end
    println("Generate instance data: ", now() - start_time)
    flush(stdout)
    
    model = Model(HiGHS.Optimizer)
    start_time = now()
    @variable(model, 0 <= x[col=1:length(edges)] <= 1.0);
    @variable(model, 0 <= w[j=1:m] <= 1.0);
    println("Create variables: ", now() - start_time)
    flush(stdout)

    start_time = now()
    @objective(model, Min, sum(c[col] .* x[col] for col = 1:length(edges)));
    println("Set objective: ", now() - start_time)
    flush(stdout)

    start_time = now()
    treatment_assignment_matrix = spzeros(n, length(edges));
    control_assignment_matrix = spzeros(m, length(edges));
    for col in 1:length(edges)
        (i,j) = edges[col]
        treatment_assignment_matrix[i, col] = 1
        control_assignment_matrix[j, col] = 1
    end

    # Matching constraints
    @constraint(model, treatment_assignment_matrix * x .== 1.0);
    @constraint(model, control_assignment_matrix * x .== w);

    # First moments are similar to treatment
    expected_first_moment = A' * ones(n) / n
    @expression(model, first_moment_residual, B' * w / n - expected_first_moment);
    @constraint(model, -epsilon .<= first_moment_residual .<= epsilon);

    # Second moments are similar to treatment
    for k = 1:d
        for l = k:d
            M_kl = sum(A[i,k] * A[i,l] for i = 1:n) / n
            @constraint(model, -epsilon <= sum(B[j,k] * B[j,l] * w[j] for j = 1:m) / n - M_kl <= epsilon)
        end
    end
    println("Create constraints: ", now() - start_time)
    flush(stdout)

    if optimize_model
        optimize!(model)
        solution_summary(model)
    end

    return model
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
            help = "The random seed used to generate the instance."
            arg_type = Int64
            default = 1
        "--epsilon"
            help = "Tolerance for covariate distributions of control to match the treatment."
            arg_type = Float64
        "--num_treatment_samples"
            arg_type = Int64
            default = 5000
        "--num_control_samples"
            arg_type = Int64
            default = 20000
        "--num_covariates"
            arg_type = Int64
            default = 8
        "--num_edges_per_treatment"
            arg_type = Int64
            default = 10
        "--control_shift_magnitude"
            arg_type = Float64
            default = 0.1
        "--optimize_model"
            help = "If this flag is set then the model will be optimized with HiGHS. 
            Only use this option is the model is small enough for HiGHS to optimize."
            action = :store_true
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
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    println("building model ...")
    flush(stdout)

    Random.seed!(parsed_args["seed"])

    @time "Build model" model = build_synthetic_design_matching_problem(
        parsed_args["epsilon"],
        parsed_args["num_treatment_samples"],
        parsed_args["num_control_samples"],
        parsed_args["num_covariates"],
        parsed_args["num_edges_per_treatment"],
        parsed_args["control_shift_magnitude"],
        parsed_args["optimize_model"]
    )
    flush(stdout)

    if parsed_args["rescale_model"]
        println("rescaling model ...")
	flush(stdout)
        @time "Rescale model" model = rescale_instance(lp_matrix_data(model))
	flush(stdout)
    end

    println("writing model to file ...")
    flush(stdout)
    @time "Write model" write_to_file(model, parsed_args["output_file"])
    flush(stdout)
end

main()
