using ArgParse
using JuMP
include("utils.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file"
        help = "This is the mps file that will be read and rescaled."
        arg_type = String
        required = true
        "--output_file"
        help = "This is the location that the mps file will be written to."
        arg_type = String
        required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg, val) in parsed_args
        println("  $arg  =>  $val")
    end

    original_model = read_from_file(parsed_args["input_file"])
    relaxed_model = relax_integrality(original_model)

    rescaled_model = rescale_instance(lp_matrix_data(relaxed_model))
    println("writing model to file ...")
    flush(stdout)
    @time "Write model" MOI.write_to_file(JuMP.backend(rescaled_model), parsed_args["output_file"])
    flush(stdout)
end

main()
