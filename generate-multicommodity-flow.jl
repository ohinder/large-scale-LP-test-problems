using JuMP
import HiGHS
using Dates
using Random
using LinearAlgebra
using SparseArrays
using Distributions
using Plots
using ArgParse
include("utils.jl")

function build_multicommodity_flow_problem(
    num_factories_per_commodity::Int64,
    num_commodities::Int64,
    num_warehouses::Int64,
    num_stores::Int64,
    additional_overtime_cost::Float64,
    folder_for_plots::String,
    optimize_model::Bool
)
    @assert num_factories_per_commodity >= 1
    @assert num_commodities >= 1
    @assert num_warehouses >= 1
    @assert num_stores >= 1
    @assert additional_overtime_cost >= 0.0

    start_time = now()
    average_demand_per_commodity = 100.0 * exp.(randn(num_commodities));
    distribution_per_commodity = Poisson.(average_demand_per_commodity);
    demand_per_commodity_store = zeros(num_commodities, num_stores);
    for s=1:num_stores
        demand_per_commodity_store[:,s] = rand.(distribution_per_commodity)
    end

    total_demand_per_commodity = demand_per_commodity_store * ones(num_stores);

    tolerance = 1e-8 # to stop possibility of numerical errors
    supply_per_factory = (1.0 + tolerance) * total_demand_per_commodity / num_factories_per_commodity;
    total_demand = sum(demand_per_commodity_store);

    warehouse_normal_capacity =  0.95 * total_demand / num_warehouses;
    additional_overtime_cost = 0.3;

    factory_locations = rand(num_commodities, num_factories_per_commodity, 2);
    warehouse_locations = rand(num_warehouses, 2);
    store_locations = rand(num_stores, 2);

    shipping_cost_from_factories_to_warehouses = zeros(num_commodities,num_factories_per_commodity,num_warehouses);
    for k=1:num_commodities
        for f=1:num_factories_per_commodity
            for w=1:num_warehouses
                shipping_cost_from_factories_to_warehouses[k,f,w] = norm(factory_locations[k,f,:] - warehouse_locations[w,:])
            end
        end
    end

    shipping_cost_from_warehouses_to_stores = zeros(num_warehouses, num_stores);
    for w=1:num_warehouses
        for s=1:num_stores
            shipping_cost_from_warehouses_to_stores[w,s] = norm(store_locations[s,:] - warehouse_locations[w,:])
        end
    end
    println("Generate instance data: ", now() - start_time)
    flush(stdout)

    if optimize_model
        model = direct_model(HiGHS.Optimizer())
    else
        model = direct_model(MOI.FileFormats.MPS.Model(generic_names = true))
    end

    start_time = now()
    @variable(model, flow_from_factories_to_warehouses[k=1:num_commodities,f=1:num_factories_per_commodity,w=1:num_warehouses] >= 0.0, set_string_name = false);
    @variable(model, flow_from_warehouses_to_stores[k=1:num_commodities,w=1:num_warehouses,s=1:num_stores] >= 0.0, set_string_name = false);
    @variable(model, warehouse_overtime_amount[j=1:num_warehouses] >= 0.0, set_string_name = false);
    println("Create variables: ", now() - start_time)
    flush(stdout)

    start_time = now()
    @objective(model, Min, 
        sum(shipping_cost_from_factories_to_warehouses[k,f,w] * flow_from_factories_to_warehouses[k,f,w]
            for w=1:num_warehouses for f=1:num_factories_per_commodity for k=1:num_commodities) +
        sum(shipping_cost_from_warehouses_to_stores[w,s] * flow_from_warehouses_to_stores[k,w,s]
            for k=1:num_commodities for w=1:num_warehouses for s=1:num_stores) +
        additional_overtime_cost * sum(warehouse_overtime_amount[w] for w=1:num_warehouses) 
    );
    println("Set objective: ", now() - start_time)
    flush(stdout)

    start_time = now()
    # flow from factories does not exceed supply
    for k in 1:num_commodities
        for f in 1:num_factories_per_commodity
            @constraint(model, sum(flow_from_factories_to_warehouses[k,f,w] 
                for w in 1:num_warehouses) <= supply_per_factory[k],
                set_string_name = false)
        end
    end

    # overtime constraint
    for w in 1:num_warehouses
        @constraint(model,
            sum(flow_from_factories_to_warehouses[k,f,w]
                for k=1:num_commodities
                for f=1:num_factories_per_commodity)
            <= warehouse_normal_capacity + warehouse_overtime_amount[w],
            set_string_name = false)
    end

    # flow into warehouses equals flow out of warehouses
    for w in 1:num_warehouses
        for k in 1:num_commodities
            @constraint(model, 
            sum(flow_from_factories_to_warehouses[k,f,w] 
                for f=1:num_factories_per_commodity) ==
                sum(flow_from_warehouses_to_stores[k,w,s] 
                    for s in 1:num_stores),
            set_string_name = false)
        end
    end

    # flow into stores meets demand
    for k=1:num_commodities
        for s=1:num_stores
            @constraint(model, 
                sum(flow_from_warehouses_to_stores[k,w,s] for w = 1:num_warehouses) ==
                    demand_per_commodity_store[k,s],
                set_string_name = false
            )
        end
    end
    println("Create constraints: ", now() - start_time)
    flush(stdout)

    if optimize_model
        println("optimizing model ...")
        optimize!(model);
        solution_summary(model)
        if folder_for_plots != ""
            mkdir(folder_for_plots)
            for k in randperm(num_commodities)[1:6]
                p = plot(ylims=(0,1.3))
                _flow_from_factories_to_warehouses = value.(flow_from_factories_to_warehouses);
                flow_scaling_factor = 20.0/maximum(_flow_from_factories_to_warehouses)
                first_line_drawn = true
                for w in 1:num_warehouses
                    for f in 1:num_factories_per_commodity
                        if _flow_from_factories_to_warehouses[k, f, w] > 1e-4 / flow_scaling_factor
                            if first_line_drawn
                                first_line_drawn = false
                                label = "amount shipped"
                            else
                                label = nothing
                            end
                            plot!(
                                [warehouse_locations[w,1], factory_locations[k,f,1]],
                                [warehouse_locations[w,2], factory_locations[k,f,2]],
                                color = :black,
                                label = label,
                                linewidth = _flow_from_factories_to_warehouses[k, f, w] * flow_scaling_factor,
                            )
                        end
                    end
                end
                
                _flow_from_warehouses_to_stores = value.(flow_from_warehouses_to_stores);
                for w in 1:num_warehouses
                    for s in 1:num_stores
                        if _flow_from_warehouses_to_stores[k, w, s] > 1e-4 / flow_scaling_factor
                            plot!(
                                [warehouse_locations[w,1], store_locations[s,1]],
                                [warehouse_locations[w,2], store_locations[s,2]],
                                color = :black,
                                label = nothing,
                                linewidth = _flow_from_warehouses_to_stores[k, w, s] * flow_scaling_factor,
                            )
                        end
                    end
                end

                scatter!(
                    factory_locations[k,:,1],
                    factory_locations[k,:,2],
                    markershape = :square,
                    markercolor = :red,
                    markersize = 6,
                    markerstrokecolor = :red,
                    markerstrokewidth = 2,
                    label = "Factory for commodity $k",
                )
                scatter!(
                    warehouse_locations[:,1],
                    warehouse_locations[:,2],
                    markershape = :circle,
                    markercolor = :blue,
                    markersize = 6,
                    markerstrokecolor = :blue,
                    markerstrokewidth = 2,
                    label = "Warehouses",
                )
                scatter!(
                    store_locations[:,1],
                    store_locations[:,2],
                    markershape = :cross,
                    markercolor = :green,
                    markersize = 3,
                    markerstrokecolor = :green,
                    markerstrokewidth = 2,
                    label = "Stores",
                )
                savefig(p, "$folder_for_plots/factory-warehouse-k=$k.pdf")
            end
        end
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
        "--num_factories_per_commodity"
            help = "The number of factories per commodity in the generated instance"
            arg_type = Int64
            default = 5
        "--num_commodities"
            help = "Number of commodities"
            arg_type = Int64
            default = 10
        "--num_warehouses"
            help = "Number of warehouses"
            arg_type = Int64
            default = 10
        "--num_stores"
            help = "Number of stores"
            arg_type = Int64
            default = 100
        "--additional_overtime_cost"
            help = "Number of stores"
            arg_type = Float64
            default = 0.3
        "--optimize_model"
            help = "If this flag is set then the model will be optimized with HiGHS. 
            Only use this option is the model is small enough for HiGHS to optimize."
            action = :store_true
        "--folder_for_plots"
            help = "This is the folder that plots of the optimal solution will be 
            placed. If this folder is not set then plots will not be created. Note that
            plots will only be created if the optimize model flag is set. Furthermore
            this folder cannot already exist as julia will create the folder."
            arg_type = String
            default = ""
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

    Random.seed!(parsed_args["seed"])

    if isdir(parsed_args["folder_for_plots"])
        throw(error("Folder $(parsed_args["folder_for_plots"]) already exists. Please choose a folder location that does not already exist."))
    end

    if parsed_args["folder_for_plots"] != "" && parsed_args["optimize_model"] == false
        throw(error("If folder_for_plots is nonempty then optimize_model flag should be set"))
    end

    @time "Build model" model = build_multicommodity_flow_problem(
        parsed_args["num_factories_per_commodity"],
        parsed_args["num_commodities"],
        parsed_args["num_warehouses"],
        parsed_args["num_stores"],
        parsed_args["additional_overtime_cost"],
        parsed_args["folder_for_plots"],
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
    @time "Write model" MOI.write_to_file(JuMP.backend(model), parsed_args["output_file"])
    flush(stdout)
end

main()
