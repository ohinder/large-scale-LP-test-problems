using Distances
import HiGHS
using LinearAlgebra
using Random
using SparseArrays
using ArgParse
using JuMP
import Logging
Logging.disable_logging(Logging.Warn)
# using Gurobi
# ENV["GRB_LICENSE_FILE"] = "/Users/haihaolu/gurobi.lic"
# gurobi_env = Gurobi.Env()

# We consider the production problem from Ben-Tal et al (2004). 

# E     -- Number of factories
# T     -- Number of time periods evenly spaced over one calendar year
# Dmin  -- Minmum demand in each period, where Dmin[t] is minimum demand
#          on period t
# Dmax  -- Maximum demand in each period, where Dmax[t] is maximum demand
#          on period t
# α     -- Production costs, where α[t][e] is production cost for
#          factory e on period t
# p     -- Maximum production per period per factory, where p[t][e] is the
#          maximum production for factory e on period t
# Q     -- Maximum production per year per factory, where Q[e] is the maximum
#          production for factory e
# Vmin  -- Minimum inventory at warehouse per period
# Vmax  -- Maximum inventory at warehouse per period
# v1    -- Initial inventory

function CreateProblemInstance(E,T,Dmin,Dmax,α,p,Q,Vmin,Vmax,v1,optimize_model)


    ###########################################################################
    # Specify the dimensions of the multi-stage robust linear optimization
    # problem.
    ###########################################################################

    # Because we have to make decisions in each period without knowledge of
    # the current demand, and because we have an offset term, there will be
    # a total of T+1 periods, in which decisions are made in periods 
    # 1,…,T and uncertainty is revealed on periods 2,…,T. The uncertainty
    # on the first period is the constant 1.

    # Store the number of constraints
    m = E + 2*E*T + 2*T

    # Store the number of decision variables in each period
    n = E

    # Store the number of uncertainty variables in each period
    # In this problem, we have the demand for each facility plus a
    # constant demand equal to one.
    d = 1

    # Store the number of time periods
    # Since we have a 0-th period, we will increment T by one
    
    # Store the number of constraints in the uncertainty set
    # In this problem, we have lower and upper constraints for the demand in
    # all of the periods
    r = 2*(T+1)


    ###########################################################################
    # Modify min and max for demands so that the first period is 1 and the rest
    # are incremented by one
    ###########################################################################

    Dmin_old = deepcopy(Dmin)
    Dmax_old = deepcopy(Dmax)
    Dmin = Vector{Float64}()
    Dmax = Vector{Float64}()
    push!(Dmin,1)
    push!(Dmax,1)
    for t=1:T
        push!(Dmin, Dmin_old[t])
        push!(Dmax, Dmax_old[t])
    end


    ###########################################################################
    # Set up optimization problem
    ###########################################################################

    model = JuMP.Model(HiGHS.Optimizer)

    # Decision variables for linear decision rule
    @variable(model, y[1:T,1:T,1:E])

    # Objective
    @variable(model, ε_pos[1:T] ≥ 0)
    @variable(model, ε_neg[1:T] ≥ 0)
    @constraint(model, Obj_Eq[s=1:T], ε_pos[s] - ε_neg[s]  == sum(sum(α[t][e]*y[t,s,e] for t=s:T) for e=1:E))
    @objective(model, Min, sum(Dmax[s]*ε_pos[s] - Dmin[s]*ε_neg[s] for s=1:T))

    # Capacity constriants
    @variable(model, γ_pos[t=1:T,s=1:t,e=1:E] ≥ 0)
    @variable(model, γ_neg[t=1:T,s=1:t,e=1:E] ≥ 0)
    @constraint(model, Cap_Eq[t=1:T,s=1:t,e=1:E], y[t,s,e] == γ_pos[t,s,e] - γ_neg[t,s,e])
    @constraint(model, Cap_UB[t=1:T,e=1:E], sum(Dmax[s]*γ_pos[t,s,e] - Dmin[s]*γ_neg[t,s,e] for s=1:t) ≤ p[t][e])
    @constraint(model, Cap_LB[t=1:T,e=1:E], sum(-Dmin[s]*γ_pos[t,s,e] + Dmax[s]*γ_neg[t,s,e] for s=1:t) ≤ 0)

    # Per-factory capacity constraints
    @variable(model, δ_pos[s=1:T,e=1:E] ≥ 0)
    @variable(model, δ_neg[s=1:T,e=1:E] ≥ 0)
    @constraint(model, Fac_Eq[s=1:T,e=1:E], sum(y[t,s,e] for t=s:T) == δ_pos[s,e] - δ_neg[s,e])
    @constraint(model, Fac_UB[e=1:E], sum(Dmax[s]*δ_pos[s,e] - Dmin[s]*δ_neg[s,e] for s=1:T) ≤ Q[e])

    # Total constraints
    @variable(model, η_pos[t=1:T,s=1:T] ≥ 0)
    @variable(model, η_neg[t=1:T,s=1:T] ≥ 0)
    @constraint(model, Tot_Eq[s=2:T,t=s:T], -1 + sum(sum( y[ℓ,s,e] for e=1:E) for ℓ=s:t) == η_pos[t,s] - η_neg[t,s])
    @constraint(model, Tot_UB[t=1:T],  v1 + sum(sum(y[s,1,e] for e=1:E) for s=1:t) + sum( Dmax[s]*η_pos[t,s] - Dmin[s]*η_neg[t,s] for s=2:t) - Dmin[t+1] ≤ Vmax)
    @constraint(model, Tot_LB[t=1:T], -v1 - sum(sum(y[s,1,e] for e=1:E) for s=1:t) + sum(-Dmin[s]*η_pos[t,s] + Dmax[s]*η_neg[t,s] for s=2:t) + Dmax[t+1] ≤ - Vmin)

    if optimize_model
        optimize!(model)
    end
    return model

end

function CreateProblemInstanceParameters(T,E,θ)
 
    # Intensity parameter
    ϕ =  [1 + 0.5*sin(2*π*(t-1) / T) for t in 1:T]

    # Uncertainty sets
    Dmin = [1000*(1-θ)*ϕ[t]/(T/24) for t=1:T]
    Dmax = [1000*(1+θ)*ϕ[t]/(T/24) for t=1:T]
   
    # Production costs
    α = [[(1.0 + (e-1) / (max(E - 1,1))) * ϕ[t]  for e in 1:E] for t in 1:T]
    #α = [[rand() for e in 1:E] for t in 1:T]

    # Maximum production per stage, per factory
    p = [[567/(T/24 * E/3) for e in 1:E] for t=1:T]  
    #p = [[567 for e in 1:E] for t=1:T]  

    # Maximum production over all stages, per factory
    Q = [13600 / (E/3) for e in 1:E] 

    # Bounds on warehouse inventory
    Vmin = 500     
    Vmax = 2000

    # Initial inventory 
    v1 = Vmin  
  
    # Return problem instance parameters
    return Dmin, Dmax, α, p, Q, Vmin, Vmax, v1
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--num_factories"
            help = "The number of factories in the generated instance"
            arg_type = Int64
            default = 10
        "--num_stages"
            help = "Number of time periods evenly spaced over one calendar year"
            arg_type = Int64
            default = 20
        "--uncertainty_level"
            help = "The uncertainty level in the model"
            arg_type = Float64
            default = 0.2
        "--output_file"
            help = "This is the location that the mps file will be written to."
            required = true
        "--optimize_model"
            help = "If this flag is set then the model will be optimized with HiGHS. 
            Only use this option is the model is small enough for HiGHS to optimize."
            action = :store_true
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
    E = parsed_args["num_factories"]
    T = parsed_args["num_stages"]
    θ = parsed_args["uncertainty_level"]
    optimize_model = parsed_args["optimize_model"]

    # Create problem instance
    Dmin, Dmax, α, p, Q, Vmin, Vmax, v1 = CreateProblemInstanceParameters(T,E,θ)
    model = CreateProblemInstance(E,T,Dmin,Dmax,α,p,Q,Vmin,Vmax,v1,optimize_model)

    # Write as a mps file
    write_to_file(model, parsed_args["output_file"])

end

main()

