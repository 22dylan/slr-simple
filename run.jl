# loading in other packages now
using Agents
using DataFrames
using GeoDataFrames
using ArchGDAL
using ProgressMeter
using CSV
using ArgParse
using Random
using Distributions
using StaticArrays    # necessary?
using StatsBase
import GeoFormatTypes as GFT

include(joinpath("model-backend", "ResidentialAgent.jl"))
include(joinpath("model-backend", "MiscFunc.jl"))
include(joinpath("model-backend", "MyModel.jl"))
include(joinpath("model-backend", "MyParcelSpace.jl"))



"""
TODOs:
    - implement HUA
    - use abmrng(model) anywhere random numbers are generated!!
    - currently only considering housing stock (HS1, HS2, HS3); note to consider more
    - revisit economic loss calculations
    - figure out invalid surge-wave archetypes
    - check on distance calculations and projections; confirm that these are accurate.
    - clean up reading agent-networks
    - try having shared replay buffer, but independent neural nets for each agent.
"""
function main(model_runname)
    parsed_args = parse_commandline()
    if parsed_args["model_runname"] != nothing
        model_runname = parsed_args["model_runname"]
    end
    println("\nRunning Model: $model_runname")

    path_to_input_dir = joinpath(pwd(), "model-runs", model_runname, "input")
    input_dict = read_dir_CSV(path_to_input_dir)                              # reading all CSV files in directory
    df_input = input_dict["input"]
    display(df_input)

    run_model(input_dict, df_input, parsed_args, path_to_input_dir, model_runname)
end


function run_model(input_dict, df_input, parsed_args, path_to_input_dir, model_runname)
    # number of simulations
    n_iterations  = read_input_file_key(df_input, "n_iterations"; dtype=Int64)
    if parsed_args["end_sim"] !== nothing
        n_iterations = parsed_args["end_sim"]
    end

    # number of years to simulate
    n_years = read_input_file_key(df_input, "n_years"; dtype=Int64)

    # random seed
    seed = read_input_file_key(df_input, "seed"; dtype=Int64)

    # slr scenario
    slr_scenario = read_input_file_key(df_input, "slr_scenario", dtype=String)

    # slr non-exceedance prob
    slr_ne = read_input_file_key(df_input, "slr_ne", dtype=String)

    # --- Running model
    for i = parsed_args["start_sim"]:n_iterations
        model = initialize(input_dict, path_to_input_dir, model_runname, seed=seed, iter=i)
        data_a, data_m, data_s = my_run!(
                                    model,
                                    n_years;
                                    adata=get_agent_save_data(),
                                    mdata=get_model_save_data(),
                                    sdata=get_space_save_data(),
                                    )

        # --- saving results for iteration
        fn_agnts = "df_agnts_$(i)_sc$(slr_scenario)_ne$(slr_ne).csv"
        write_out(data_a, model, fn_agnts)

        fn_model = "df_model_$(i)_sc$(slr_scenario)_ne$(slr_ne).csv"
        write_out(data_m, model, fn_model)

        fn_space = "df_space_$(i)_sc$(slr_scenario)_ne$(slr_ne).csv"
        write_out(data_s, model, fn_space)

        close_model!(model, model_runname)
    end
end





"""
    my_run!()
Custom model run function. 
Started from example on Agents.jl docs
Needed to customize to collect space (pacel) data during model time steps
"""
function my_run!(model,
                n;
                mdata = nothing,
                adata = nothing,
                sdata = nothing,
                )

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    df_space = init_space_dataframe(model, sdata)
    s = 0    
    while Agents.until(s, n, model)
      collect_agent_data!(df_agent, model, adata)
      collect_model_data!(df_model, model, mdata)
      collect_space_data!(df_space, model, sdata, s)
      step!(model, 1)
      s += 1
    end
    return df_agent, df_model, df_space
end


"""
    init_space_dataframe()
Function to initialize space dataframe.
Used to store space (parcel) results for output
"""
function init_space_dataframe(model::ABM, properties::AbstractArray)
    std_headers = 2

    headers = Vector{String}(undef, std_headers + length(properties))
    headers[1] = "step"
    headers[2] = "guid"

    for i in 1:length(properties)
        headers[i+std_headers] = dataname(properties[i])
    end
    types = Vector{Vector}(undef, std_headers + length(properties))
    types[1] = Int[]
    types[2] = String[]
    for (i, field) in enumerate(properties)
        types[i+2] = eltype(model.prcl_df[!,field])[]
    end

    push!(properties, :guid)
    push!(properties, :step)
    return DataFrame(types, headers)
end


"""
    collect_space_data!()
Function to collect space (parcel) data from model at each time step
"""
function collect_space_data!(df, model, properties::Vector, step::Int = 0)
    dd = model.prcl_df
    dd[!,:step] = fill(step, size(dd,1))
    dd = dd[!,properties]
    append!(df, dd)
    return df
end


"""
    parse_commandline()
Function to pass option values into the model from the command line
Created to start simulation from user-specified iteration
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--start_sim"
            help = "option to start simulation at iteration other than 1"
            default = 1
            arg_type = Int
            required = false
        "--end_sim"
            help = "option to end simulation at iteration other than 1"
            default = nothing
            arg_type = Int
            required = false
        "--model_runname"
            help = "option to run specific model runname"
            default = nothing
            arg_type = String
            required = false
    end
    return parse_args(s)
end


# ------------
model_runnames = [
                "Status-quo-2",             # status quo - this runs for all of Galveston (21,353 buildings)
                # "testbed",                  # testbed    - this runs for a subset of Galveston (1,284 buildings)
                ]


for model_runname in model_runnames
    main(model_runname)
end