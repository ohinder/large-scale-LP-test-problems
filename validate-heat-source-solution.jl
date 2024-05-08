"""
This file takes in the ground truth heat solution to temperature profile to the heat source problem
and the solution produced by PDLP and validates that they are close. This can demonstrate that PDLP 
has succesfully recovered the true solution.
"""

import HDF5
using JuMP, LinearAlgebra, ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ground_truth_file"
        help = "This is the file that includes the data for the ground truth value of the temperature profile u.
        This is a HDF5 file produced when the instance is generated (using generate-heat-source-location.jl)."
        arg_type = String
        required = true
        "--solution_file"
        help = "This is the location that solution from the solver was written. It accepts .sol files 
        outputted by PDLP."
        required = true
        arg_type = String
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg, val) in parsed_args
        println("  $arg  =>  $val")
    end

    ground_truth_file = parsed_args["ground_truth_file"] #"/Users/Oliver/Downloads/easy/temperature_ground_truth1-names.txt"
    u_true = HDF5.h5read(ground_truth_file, "u_true");

    solution_file = parsed_args["solution_file"] #"/Users/Oliver/Downloads/easy/heat-source-instance-easy-names.sol"

    file = open(solution_file)
    total_variables = countlines(file) - 1
    seekstart(file)
    solution_vector = Array{Float64,1}(undef, total_variables)
    i = 0
    for line in eachline(file)
        if i > 0
            varname, value = split(line, " ")
            solution_vector[i] = parse(Float64, value)
        end
        i += 1
    end
    close(file)

    # solve equation 
    # n^3 + (n-2)^3 = total_variables to compute the number of heat variables 
    n = Int(ceil((total_variables / 2)^(1/3)))
    @assert n^3 + (n-2)^3 == total_variables

    u_solver = solution_vector[1:n^3]
    @assert length(u_solver) == length(u_true[:])

    println("Relative solution error:")
    @show norm(u_solver - u_true[:]) / norm(u_true[:])
end

main()