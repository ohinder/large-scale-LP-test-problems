using JuMP
import HiGHS
using ArgParse
include("utils.jl")
# References:
# - https://www.opt.math.tugraz.at/~cela/papers/qap_bericht.pdf
# - https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/
# - P. Adams and T. A. Johnson, Improved linear programming-based lower bounds for the quadratic assignment problem, in Quadratic Assignment and Related Problems, P. M. Pardalos and H. Wolkowicz, eds., DIMACS Series on Discrete Mathematics and Theoretical Computer Science 16, 1994, 43–75, AMS, Providence, RI.

function read_QAP_library_problem(filename::String)
    # Assumes file is in the format:
    # problem size  
    #
    # Flow matrix 
    #
    # Distance matrix

    # Note it expects a new line between problem size, flow matrix and distance matrix 
    # Some of the QAPLIB problems have this space and some don't -- you may have to manually 
    # add the space.

    # chatGPT helped write this function

    # Open the file
    open(filename, "r") do file
        # Read the number n
        n = parse(Int, readline(file))
        
        # Skip the empty line after n
        readline(file)
                
        B = zeros(Int, n, n) # distance matrix
        
        # Read matrix A
        A_vec = Int[] 
        while length(A_vec) < n^2 
            append!(A_vec, parse.(Int, split(readline(file))))
        end 
        A = reshape(A_vec, (n,n))

        # Skip the empty line after matrix A
        readline(file)
        
        # Read matrix B
        B_vec = Int[] 
        while length(B_vec) < n^2 
            append!(B_vec, parse.(Int, split(readline(file))))
        end 
        B = reshape(B_vec, (n,n))

        return A, B
    end
end

function write_qap_problem_in_terri_johnsons_format(A::Matrix, B::Matrix, filename::String)
    # Only works for problems with symmetric flow and distance matricies 
    @assert size(A, 1) == size(A, 2)
    @assert size(B, 1) == size(B, 2)
    @assert size(A, 1) == size(B, 1)

    n = size(A,1)
    matrix = zeros(Int, n, n)
    for i = 1:n 
        @assert A[i,i] == 0
        @assert B[i,i] == 0
        for j = 1:n 
            if i > j
                matrix[i,j] = A[i,j] # flow in lower half
                @assert A[i,j] == A[j,i]
            elseif i < j
                matrix[i,j] = B[i,j] # distances in upper half
                @assert B[i,j] == B[j,i]
            end 
        end
    end 

    # written with chatGPT
    open(filename, "w") do file
        # Write the matrix size
        println(file, n)
        
        # Write the matrix elements
        for row in 1:n
            for col in 1:n
                print(file, matrix[row, col])
                if col < n
                    print(file, ",")
                end
            end
            println(file)
        end
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file"
            help = "The input file from QAPLIB. You can download the QAPLIB files 
            from https://coral.ise.lehigh.edu/data-sets/qaplib/. Note: it should 
            be a symmetric instance and expects the input files to have a blank line before
            each matrix (not all the data files have this but you can add it in manually)."
            required = true
        "--output_file"
            help = "This is the location that the file in TJ format will be written to."
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

    A, B = read_QAP_library_problem(parsed_args["input_file"])     
    write_qap_problem_in_terri_johnsons_format(A, B, parsed_args["output_file"])
end

main()
