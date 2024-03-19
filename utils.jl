using JuMP, LinearAlgebra
# Returns a new rescaled jump model.
#  If this is true we rescale the model. This rescaling makes 
# the objective coefficients all 0, 1 or -1. The rescaling makes all the right hand 
# sides either -1, 0 or 1, except in the case that a constraint has two nonzero finite, 
# nonmatching lower and upper bounds. In this case we rescale such that the maximum 
# absolute value of the right hand side is one.
# Note: variable and constraint name information is currently thrown away in this process
function rescale_instance(data::JuMP.LPMatrixData)
    rescaling_rhs = ones(length(data.b_lower))
    for i = 1:length(data.b_lower)
        if abs(data.b_lower[i]) == Inf && 0.0 < abs(data.b_upper[i]) < Inf
            rescaling_rhs[i] = 1.0 / abs(data.b_upper[i])
        elseif abs(data.b_upper[i]) == Inf && 0.0 < abs(data.b_lower[i]) < Inf
            rescaling_rhs[i] = 1.0 / abs(data.b_lower[i])
        elseif 0.0 < max(abs(data.b_lower[i]), abs(data.b_upper[i])) < Inf
            rescaling_rhs[i] = 1.0 / max(abs(data.b_lower[i]), abs(data.b_upper[i]))
        end
    end

    rescaling_obj = ones(length(data.c))
    for i = 1:length(data.c)
        if 0.0 < abs(data.c[i]) < Inf
            rescaling_obj[i] = 1.0 / abs(data.c[i])
        end
    end

    # perform rescaling on data
    b_lower = data.b_lower .* rescaling_rhs
    b_upper = data.b_upper .* rescaling_rhs
    
    x_lower = data.x_lower ./ rescaling_obj
    x_upper = data.x_upper ./ rescaling_obj
    c = data.c .* rescaling_obj

    A = (Diagonal(rescaling_rhs) * data.A) * Diagonal(rescaling_obj)

    # build jump model with rescaled data
    rescaled_model = Model()
    @variable(rescaled_model, x_lower[i] .<= x[i=1:length(x_lower)] .<= x_upper[i], set_string_name=false)
    @objective(rescaled_model, data.sense, sum(c[i] * x[i] for i=1:length(x_lower)))
    @constraint(rescaled_model, b_lower .<= A * x .<= b_upper, set_string_name=false)
    return rescaled_model
end

# help description for the rescale_model flag.
# this is here because it is used in several different models
RESCALE_MODEL_HELP_DESCRIPTION = "If this is true we rescale the model. This rescaling makes 
            the objective coefficients all 0, 1 or -1. The rescaling makes all the right hand 
            sides either -1, 0 or 1, except in the case that a constraint has two nonzero finite, nonmatching lower 
            and upper bounds. In this case we rescale such that the maximum absolute value of the right hand side is 
            one. The purpose of this change is to ensure that the termination criteria is consistent across 
            solvers. Modifying this option will thus affect the difficulty of problem instances. 
            Also, this option may create many small nonzero coefficients in the A matrix. These should not
            be rounded to zero even if they are less than say 1e-8 because it may cause inaccurate 
            solutions to problems with many constraints or variables."
