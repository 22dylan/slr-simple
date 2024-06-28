Base.@kwdef mutable struct Parameters
    # model values
    input_df::DataFrame
    input_dir::String
    output_dir::String
    tick::Int64 = 1

    iteration::Int64 = 1
    n_iterations::Int64

    n_years::Int64
    slr_scenario::String
    slr_exposure::DataFrame
    slr_exit_tt::DataFrame
    slr_utmb_tt::DataFrame
    slr_elec_acc::DataFrame
    year::Int64
    start_year::Int64
    progress_bar::Progress
    n_prcls::Int64
    prcl_df::DataFrame
    pos_to_guid_idx::Dict
    agent_view_radius::Int64=11
    n_unoccupied::Int64=0
    n_occupied::Int64=0

end


"""
    initialize(input_dict, path_to_input_dir; seed, iter)
Sets up the initial model using the input directory
Sets up the parcel space
Assigns agents to the model
"""
function initialize(input_dict, path_to_input_dir, model_runname; seed=1337, iter=1)
    seed = seed + iter 
    rng = Random.MersenneTwister(seed)   # setting seed

    # preparing parcel dataframe and converting to julia
    cell_size = 25  # cell size in meters
    prcl_df = prepare_parcel_df(path_to_input_dir, seed=seed, cell_size=cell_size, rotate=true)

    # setting up ABM space as grid
    gridsize = (div(47_000,cell_size), div(7_000,cell_size))
    space = GridSpace(gridsize; periodic=false)       # galveston is approx. 47,000m by 7,000m

    # getting model property dictionary
    properties = setup_model_properties(model_runname, input_dict, prcl_df, iter, path_to_input_dir, rng)

    # setting up the model
    model = ABM(
                ResidentialAgent,
                space;
                rng=rng,
                properties=properties,
                model_step! =model_step!,  # note: a space is needed after "model_step!", otherwise this is interpreted as "model_step !="
                )

    # --- adding agents to model; either in parcel or general environment
    AddAgentsToModel!(model, prcl_df)
    UpdateModelCounts!(model)
    return model
end




"""
    setup_model_properties(input_dict::Dict, prcl_df::DataFrame, iter::Int64)
setting up model properties that go into the ABM.
input_dict: input dictionary; contains dataframes of CSV files
prcl_df: parcel dataframe that was previously prepared
iter: iteration number in outer ABM cycle
"""
function setup_model_properties(model_runname::String, input_dict::Dict, prcl_df::DataFrame, iter::Int64, path_to_input_dir::String, rng::AbstractRNG)
    input_df = input_dict["input"]

    n_iterations = read_input_file_key(input_df, "n_iterations"; dtype=Int64)
    n_years = read_input_file_key(input_df, "n_years"; dtype=Int64)
    slr_scenario = read_input_file_key(input_df, "slr_scenario", dtype=String)
    model_start_year = read_input_file_key(input_df, "model_start_year")

    progress_bar = ProgressBar(iter, n_iterations, n_years)
    update!(progress_bar)

    slr_exposure, slr_exit_tt, slr_utmb_tt, slr_elec_acc = read_slr_exposure(input_df, order=prcl_df.guid, rng=rng)
    pos_to_guid_idx = setup_pos2guididx(prcl_df)

    path_to_output_dir = joinpath(pwd(), "model-runs", model_runname, "output")
    makedir(path_to_output_dir)

    # ----------------------------------------
    properties = Parameters(
        input_df=input_df,
        input_dir=path_to_input_dir,
        output_dir=path_to_output_dir,
        n_iterations=n_iterations,
        n_years=n_years,
        slr_scenario=slr_scenario,
        slr_exposure=slr_exposure,
        slr_exit_tt=slr_exit_tt,
        slr_utmb_tt=slr_utmb_tt,
        slr_elec_acc=slr_elec_acc,
        year=model_start_year,
        start_year=model_start_year,
        progress_bar=progress_bar,
        n_prcls=size(prcl_df, 1),
        prcl_df=prcl_df,
        pos_to_guid_idx=pos_to_guid_idx,
        )
    return properties
end


"""
    setup_pos2guididx(prcl_df::DataFrame)
setting up a dictionary to map pos (tuple) to row in prcl_df
"""
function setup_pos2guididx(prcl_df::DataFrame)
    pos_to_guid_idx = Dict()
    for pos in unique(prcl_df.pos)
        p_ = prcl_df[findall(==(pos),prcl_df.pos),:]
        pos_to_guid_idx[pos] = p_.row
    end
    return pos_to_guid_idx
end



"""
    AddAgentsToModel!(model::ABM, parcel_df::DataFrame)
adds agents to the model to initialize
"""
function AddAgentsToModel!(model::ABM, parcel_df::DataFrame)
    add_agents!(model, parcel_df)
    setup_initial_states!(model)
end


"""
function add_agents!(model, parcel_df)
"""
function add_agents!(model::ABM, parcel_df::DataFrame)
    for i = 1:size(parcel_df,1)
        id = i
        i !=1 && (id=next_avail_id(model))   # if id!=1, then get the next available id
        p = parcel_df[i,:]                              # parcel information
        s = define_state(ResidentialAgent)
        r = div(size(s,1)-1,2)
        alphas = alpha_calc(model,3)
        agent = ResidentialAgent(
                            id,                         # agent id (required)
                            p.pos,                      # agent pos as tuple (required)
                            p.guid,                     # agent's current  position (guid)
                            p.row,                      # pos_idx; position index - used to quickly look up from prcl_df
                            :nothing,                   # action
                            r,                          # agent view radius
                            s,                          # agent state; captures neighbors and neighborhood exposure
                            zeros(MArray{Tuple{11,11},Float32}), # original neighborhood
                            alphas
                            )
        add_agent!(agent, agent.pos, model)
    end
end

"""
    setup_initial_states!(model)
gets each agent's initial state
updates agent.state, counts number of neighbors at model start;
this done is after adding all agents to the model; that way agents see who's around them
"""
function setup_initial_states!(model)
    for id in random_ids(model)
        model[id].neighborhood_original = get_neighborhood(model[id], model)    # get each agent's original neighborhood before any migration occurs
        update_state!(model[id], model)                                         # setting up the agent's initial state
    end
end


"""
    read_slr_exposure(input_df)
reads slr exposure data for the slr scenario in input.csv
The slr exposure data was pre-computed and is a CSV file containing each guid (rows)
and year of exposure (columns)
"""
function read_slr_exposure(input_df; order::Vector{String}, rng::AbstractRNG)
    sc = read_input_file_key(input_df, "slr_scenario", dtype=String)
    ne = read_input_file_key(input_df, "slr_ne", dtype=String)
    path_to_slr_scns = joinpath(pwd(), "slr-scenarios")
    path_to_exposure = joinpath(pwd(), "slr-scenarios", "building-exposure", "nTimesExp_years_sc$(sc)_ne$(ne).csv")
    path_to_exit_tt = joinpath(path_to_slr_scns, "galveston-exit", "nTTIncrease_years_sc$(sc)_ne$(ne).csv")
    path_to_utmb_tt = joinpath(path_to_slr_scns, "utmb-hospital", "nTTIncrease_years_sc$(sc)_ne$(ne).csv")
    path_to_elec_acc = joinpath(path_to_slr_scns, "electricity-access", "nNoAccess_years_sc$(sc)_ne$(ne).csv")

    slr_exposure = read_csv(path_to_exposure)
    slr_exit_tt = read_csv(path_to_exit_tt)
    slr_utmb_tt = read_csv(path_to_utmb_tt)
    slr_elec_acc = read_csv(path_to_elec_acc)

    slr_exposure = slr_exposure[indexin(order, slr_exposure.guid),:]
    slr_exit_tt  = slr_exit_tt[indexin(order, slr_exit_tt.guid),:]
    slr_utmb_tt  = slr_utmb_tt[indexin(order, slr_utmb_tt.guid),:]
    slr_elec_acc  = slr_elec_acc[indexin(order, slr_elec_acc.guid),:]

    # converting each column from Float64 to Int64
    convert_cols_int!(slr_exposure)
    convert_cols_int!(slr_exit_tt)
    convert_cols_int!(slr_utmb_tt)
    convert_cols_int!(slr_elec_acc)
    return slr_exposure, slr_exit_tt, slr_utmb_tt, slr_elec_acc
end


function convert_cols_int!(df)
    for col in names(df)
        if col == "guid"
            continue
        end
        df[!,col] = convert.(Int64,df[!,col])
    end
end


"""
    model_step!(model)
function for custom model step
checks time in model.
"""
function model_step!(model::ABM)
    if model.tick < model.n_years         # pre-hazard
        AllAgentsStep!(model)
    else
        AllAgentsCloseStep!(model)
    end
    UpdateModelCounts!(model)
    model.tick += 1
    model.year += 1
    next!(model.progress_bar)
end


function AllAgentsStep!(model::ABM)
    # loop through agents to get/perform action
    for id in random_ids(model)
        agent_step!(model[id], model)
    end
    return model
end

function random_ids(model::ABM)
    ids = collect(allids(model))         # getting ids of agents
    ids = shuffle(abmrng(model), ids)    # shuffling ids
    return ids
end


"""
    UpdateModelCounts!(model::ABM)
updating counts in the model.
Inlcudes number of unoccupied and occupied parcels.
TODO: currently considering parcels with csv_sector as HS1, HS2, HS3; need to consider others
"""
function UpdateModelCounts!(model::ABM)
    model.n_unoccupied = count_unoccupied(model)
    model.n_occupied = count_occupied(model)
end

"""
    count_unoccupied(model::ABM)
    count_occupied(model::ABM)
function to count number of occupied/unoccupied parcels
"""
count_occupied(model::ABM) = sum(Int64, model.prcl_df.occupied)
count_unoccupied(model::ABM) = sum(Int64, 1 .- model.prcl_df.occupied)



"""
    close_model!(model::ABM)
closes up the model; placeholder function for now
perform final agent steps and model counts
"""
function close_model!(model::ABM, model_runname::String)
    return
end

function AllAgentsCloseStep!(model::ABM)
    ids = collect(allids(model))
    ids = shuffle(abmrng(model), ids)
    for id in ids
        agent_close_step!(model[id], model)
    end
end

# function merge_dmg_with_prcl_df!(model::ABM)
#     dmg_results = read_csv(joinpath("temp", "BldgDmg-Wind.csv"))
#     model.prcl_df = innerjoin(model.prcl_df, dmg_results[:, [:guid, :DS_0, :DS_1, :DS_2, :DS_3]], on=:guid)
#     rename!(model.prcl_df, Dict(:DS_0=>:DS_0_Wind, :DS_1=>:DS_1_Wind, :DS_2=>:DS_2_Wind, :DS_3=>:DS_3_Wind))

#     dmg_results = read_csv(joinpath("temp", "BldgDmg-Flood.csv"))
#     model.prcl_df = innerjoin(model.prcl_df, dmg_results[:, [:guid, :DS_0, :DS_1, :DS_2, :DS_3]], on=:guid)
#     rename!(model.prcl_df, Dict(:DS_0=>:DS_0_Flood, :DS_1=>:DS_1_Flood, :DS_2=>:DS_2_Flood, :DS_3=>:DS_3_Flood))

#     dmg_results = read_csv(joinpath("temp", "BldgDmg-SurgeWave.csv"))
#     model.prcl_df = innerjoin(model.prcl_df, dmg_results[:, [:guid, :DS_0, :DS_1, :DS_2, :DS_3]], on=:guid)
#     rename!(model.prcl_df, Dict(:DS_0=>:DS_0_SurgeWave, :DS_1=>:DS_1_SurgeWave, :DS_2=>:DS_2_SurgeWave, :DS_3=>:DS_3_SurgeWave))
# end

function read_input_file_key(input_df::DataFrame, key::String; dtype::DataType=Int64)
    v = input_df[input_df[:,"Variable"] .== key, "Value"][1]
    (dtype==Bool) && (return parse_bool(v))

    if dtype == String
        return v
    else
        v = parse(dtype, v)
    end
    return v
end

function parse_bool(v)
    v = parse(Int64, v)
    v = Bool(v)
    return v
end

"""
    ProgressBar(i, iters, n_years)
prints status of model to terminal
"""
function ProgressBar(i, iters, n_years; color=:cyan, description="Iteration")
    p = Progress(n_years, 
            desc="$(description): $(i)/$(iters)",
            # desc="Iteration: $(i)/$(iters) | ",
            barlen=30, 
            color=color
        )
    return p
end