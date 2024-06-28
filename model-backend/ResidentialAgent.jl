const State = MArray{Tuple{11,11,2},Float32}    # state; using static array (11x11 grids by 2 metrics (number neighbors that have left, neighbor exposure))
@agent struct ResidentialAgent(GridAgent{2})
    # id::Int64                 # agent id; requied by gridspace; automatically included in @agent
    # pos::Tuple(Int64,Int64)   # agent position; required by gridspace; automatically included in @agent
    pos_guid::String            # guid that agent is associated with or "none"
    pos_idx::Int                # agent position as index; used for quicker lookups
    action::Symbol              # action chosen
    view_radius::Int64          # radius, r, that agent can view; 
    state::State                # agent state defined as mutable static array
    neighborhood_original::MArray{Tuple{11,11}} # number of neighbors at start of simulation'
    alphas::Vector{Float64}   # alpha values

end

define_state(::Type{ResidentialAgent}) = State(zeros(Float32, size(State)))

function update_state!(agent::ResidentialAgent, model::ABM; dt::Int64=0)
    update_state_migration!(agent, model)
    update_state_exposure!(agent, model, dt)
end

function update_state_migration!(agent::ResidentialAgent, model::ABM)
    s = get_neighborhood(agent, model)
    agent.state[:,:,1] = agent.neighborhood_original - s
end

function get_neighborhood(agent::ResidentialAgent, model::ABM)
    state = zeros(Float32, (model.agent_view_radius, model.agent_view_radius))  # pre-allocating space
    center_pos = (agent.view_radius+1, agent.view_radius+1)                     # center position in state matrix
    nearby_pos = nearby_positions(agent.pos, model, agent.view_radius)          # getting nearby positions in ABM model
    for near_pos in nearby_pos                                                  # loop through nearby positions
        d_pos = near_pos .- agent.pos                                           # difference between nearby position and current position
        idx = (center_pos[2]-d_pos[2], center_pos[1]+d_pos[1])                  # idx in state is (row,col); confusing because rows are changes in "y-direction", cols are changes in "x-direction"
        if haskey(model.pos_to_guid_idx, near_pos)                              # if position has a parcel in it
            id_pos = ids_in_position(near_pos, model) 
            (id_pos == Int64[]) && (state[idx...]=0f0)                          # if no agent in position, then 0
            (id_pos != Int64[]) && (state[idx...]=convert(Float32, length(id_pos)))  # if agent in position, then length of ids in pos
        else
            state[idx...] = 0f0                                                 # if no parcel in position, then 0.0
        end
    end
    # updating state of position itself
    state[center_pos...] = convert(Float32, length(ids_in_position(agent.pos, model)))
    return state
end

function update_state_exposure!(agent::ResidentialAgent, model::ABM, dt::Int64)
    state = zeros(Float32, (model.agent_view_radius, model.agent_view_radius))  # pre-allocating space
    center_pos = (agent.view_radius+1, agent.view_radius+1)                     # center position in state matrix
    nearby_pos = nearby_positions(agent.pos, model, agent.view_radius)          # getting nearby positions in ABM model
    year_string = "_$(model.year+dt)"
    for near_pos in nearby_pos                                                  # loop through nearby positions
        d_pos = near_pos .- agent.pos                                           # difference between nearby position and current position
        idx = (center_pos[2]-d_pos[2], center_pos[1]+d_pos[1])                  # idx in state is (row,col); confusing because rows are changes in "y-direction", cols are changes in "x-direction"
        if haskey(model.pos_to_guid_idx, near_pos)                              # if position has a parcel in it
            e = 0f0                                                             # pre-allocating exposure; 0f0 is a 0.0 in Float32
            prcl_indices = model.pos_to_guid_idx[near_pos]                      # get prcl_df indices of parcels in cell
            for prcl_i in prcl_indices                                          # loop through parcels in the cell
                e_ = get_prcl_info(model.slr_exposure, prcl_i, year_string)     # getting slr exposure of parcel
                (e_>e) && (e=e_)                                                # taking max exposure of parcel in cell
            end
            e = e/365                                                           # normalzing to percent of year
            (e>1) && (e=1)                                                      # if greater than 1, set to 1
            state[idx...] = convert(Float32,e)
        else
            state[idx...] = 0f0                                                 # if no parcel in position, then 0.0
        end
    end
    # updating state of position itself
    e = get_prcl_info(model.slr_exposure, agent.pos_idx, year_string)
    e = e/365
    (e>1) && (e=1)
    state[center_pos...] = convert(Float32, e)
    agent.state[:,:,2] = state
end

function get_exposure(pos::Tuple, model::ABM, tick::Int64, year::Int64; norm::Bool=true)
    # testing out letting the agent have memory of parcel that they left
    pos_idx = model.pos_to_guid_idx[pos]
    year_prior_string = "none"
    year_string = "_$(year)"
    (tick>1) && (year_prior_string = "_$(year-1)")

    days_exp = get_prcl_info(model.slr_exposure, pos_idx[1], year_string)
    if year_prior_string != "none"
        days_exp_prior = get_prcl_info(model.slr_exposure, pos_idx[1], year_prior_string)
        days_low_access_exit_prior = get_prcl_info(model.slr_exit_tt, pos_idx[1], year_prior_string)
    else
        days_exp_prior = 0.0
        days_low_access_exit_prior = 0.0
    end

    days_low_access_exit = get_prcl_info(model.slr_exit_tt, pos_idx[1], year_string)
    days_low_access_utmb = get_prcl_info(model.slr_utmb_tt, pos_idx[1], year_string)

    (days_exp>365) && (days_exp=365)                            # if 366, assuming 365
    (days_exp_prior>365) && (days_exp_prior=365)                # if 366, assuming 365
    (days_low_access_exit>365) && (days_low_access_exit=365)    # if 366, assuming 365
    (days_low_access_utmb>365) && (days_low_access_utmb=365)    # if 366, assuming 365

    (norm==true) && (days_exp=days_exp/365)             # option to normalize
    (norm==true) && (days_exp_prior=days_exp_prior/365) # option to normalize
    (norm==true) && (days_low_access_exit=days_low_access_exit/365) # option to normalize
    (norm==true) && (days_low_access_exit_prior=days_low_access_exit_prior/365) # option to normalize
    (norm==true) && (days_low_access_utmb=days_low_access_utmb/365) # option to normalize

    return days_exp, days_exp_prior, days_low_access_exit, days_low_access_exit_prior, days_low_access_utmb
end

"""
    act!(agent, model, action)
residential agent action,
if agent has already left, then action is nothing.
  note that agent still collects reward if they've left
"""
function act!(agent::ResidentialAgent, model::ABM, action::Symbol)
    (action == :nothing) && (return)
    (action == :leave) && (agent_leaves!(agent, model))
end

"""
    agent_leaves!(agent::ResidentialAgent, model::ABM)
function for removing agent from model
updates the parcel dataframe such that occupied column is 0.0
removes the agent from the ABM
"""
function agent_leaves!(agent::ResidentialAgent, model::ABM)
    model.prcl_df[agent.pos_idx[1], :occupied] = 0.0    # updating parcel dataframe so that :occupied is 0.0
    remove_agent!(agent, model)
end

"""
    agent_close_step!(agent::ResidentialAgent, model::ABM)
place-holder function for agent's final step; not needed now, but could be later
"""
function agent_close_step!(agent::ResidentialAgent, model::ABM)
    nothing
end

################################################################################

"""
    agent_step!(agent::ResidentialAgent, model::ABM)

Nat: you'll likely be starting from this function and working in this file

This function is repeatedly called for all ResidentialAgents and defines how the 
  agents behave.
Each agent has a state (`agent.state`) and can take an action (agent.action) at 
  each time step.
The agents will be performing actions based on their state
There are 75 time steps (corresponding to 2025-2100) in one iteration. 
We usually run multiple iterations.

## `agent.state`
- Each agent's state is represented as two matrices (stored in an 11x11x2 array).
- Each of these two matrices represents different attributes of that agent's 
  neighborhood. 
- The agent's own location is at the center of the matrix.

- The first 11x11 matrix (`agent.state[:,:,1]`) defines migration in the 
  agent's neighborhood. for example, if `agent.state[1,3,1] = 1.0`, this means 
  that another agent in that location has migrated b/c of sea level rise (SLR).
- An example (5x5) matrix, could look like this; meaning that 4 neighbors have 
  left
        [0, 0, 0, 0, 0
         0, 1, 0, 0, 0
         0, 1, 0, 1, 0
         0, 0, 0, 0, 0
         0, 0, 0, 0, 1]

- The second 11x11 matrix (`agent.state[:,:,2]`) defines how often other buildings 
  in the neighborhood are exposed to SLR. 
- This is represented as a percentage of the year (e.g., 0.3 -> that building is 
  exposed 30% of the year).
- An example (5x5) matrix, could look like this; here, there are two buildngs that
  are exposed 30% and 10% of the year
        [0, 0,   0, 0, 0
         0, 0.3, 0, 0, 0
         0, 0.1, 0, 0, 0
         0, 0,   0, 0, 0
         0, 0,   0, 0, 0]

## `agent.action`
- This is used to define the agent's action
- We'll start with considering two actions: (1) do nothing, and (2) migrate
    - if the agent is doing nothing, this is represented as `agent.action = :nothing`
    - if the agent decides to migrate, this is changed to `agent.action = :leave`

- You'll be playing around with how we define when agent's perform an action
- There's a function called `act!()` that will tell the rest of the model what to 
  do based on `agent.action`

Below are some thoughts on first steps in this version of the ABM:
    - create simple rules where agents make decisions only on what's happening to
      them. 
        - the agents ignore their neighbors here. 
        - an example rule could be if the agent is exposed >30% of the year, then
          they migrate/leave
        - this is a great starting point since it's very simple. 
    - simple rules where agents make decisions based on what's happening to both 
      them and their neighbors. 
        - this is a next step from the above and considers all of `agent.state`
        - an example rule could be the agent decides to leave if:
            - they are exposed >30% of the year, OR
            - 30% of their neighborhood has left

Below are some thoughts on next steps after the above two are done:
    - have agents make decisions using utility funcitons 
        - https://en.wikipedia.org/wiki/Utility
        - Utility functions are commonly used in economics (e.g., Cobb-Douglas)
        - If we do this, I'll need to think a little more about how exactly this 
          should be written/coded; 
        - We can also talk to Jenn Helgeson, an economist at NIST.

    - add in how often the agent loses electricity and increases in travel time.
        - I can do the adding in part to the model since I have these results 
        - We can talk about how you can use it for agent decisions.

    - add distance-decay function for agent neighborhood
        - https://en.wikipedia.org/wiki/Distance_decay
        - this would make it so that agent's don't care as much about what's 
          happening greater distances from them in their neighborhood
"""

function test_one(agent::ResidentialAgent, model::ABM)
    update_state!(agent, model)         # updates agent state for this time step
    # Extract the exposure of the agent's own position from the state
    center_pos = (agent.view_radius + 1, agent.view_radius + 1)
    exposure = agent.state[center_pos..., 2]

    # Define the thresholds for leaving
    threshold_1 = 0.1
    threshold_2 = 0.3
    threshold_3 = 0.5

    # Determine action based on exposure

    
    if exposure > threshold_3
        agent.action = :leave
    else
        agent.action = :nothing
    end
end

function test_two(agent::ResidentialAgent, model::ABM)
    # Constants
    PR = 0.1
    PRTWO = 0.05
    PRThree = 1.0 

    # Update the agent's state for the current time step
    update_state!(agent, model)
    
    # Extract the exposure of the agent's own position from the state
    center_pos = (agent.view_radius + 1, agent.view_radius + 1)
    ER =-1.0* (agent.state[center_pos..., 2])
    
    # Calculate the number of neighbors at the start of the simulation
    num_neighbors_start = sum(agent.neighborhood_original)

    # Calculate the number of neighbors at the current time step
    num_neighbors_current = sum(agent.state[:, :, 1])

    # Calculate NR using the given formula
    NR = 0.05 * (1 - (num_neighbors_current / num_neighbors_start))

    # Calculate the reward value (r)
    r = PRThree + NR + ER

    # Determine action based on the value of r
    if r <= 0
        agent.action = :leave
    else
        agent.action = :nothing
    end
end


function utility_theory(agent::ResidentialAgent, model::ABM)
    # Calculate UStay for the current agent
    UStay_self = 100 ^ agent.alphas[1]

    # Initialize 2nd and 3rd alpha values
    alpha_2 = agent.alphas[2]
    alpha_3 = agent.alphas[3]

    # Calculate the number of neighbors at the start of the simulation
    num_neighbors_start = sum(agent.neighborhood_original)

    # Calculate the number of neighbors at the current time step
    num_neighbors_current = sum(agent.state[:, :, 1])

    # Calculate P_neighbor using the given formula
    P_neighbor = 100 * (num_neighbors_current / num_neighbors_start)

    # Calculate Ustay_neighbor
    UStay_neighbor = P_neighbor ^ alpha_2

    # Determine UStay
    UStay = UStay_self * UStay_neighbor

    # Update state
    update_state!(agent, model)  # Ensure this function is defined elsewhere in your code

    # Extract the exposure of the agent's own position from the state
    center_pos = (agent.view_radius + 1, agent.view_radius + 1)
    exposure = agent.state[center_pos..., 2] * 100

    # Calculate ULeave
    ULeave = exposure ^ alpha_3

    # Agent decision logic
    if UStay > ULeave
        agent.action = :nothing
    elseif ULeave > UStay
        agent.action = :leave
    end

    return agent
end

function agent_step!(agent::ResidentialAgent, model::ABM)
   
   
    # test_two(agent,model)
    # test_one(agent,model)
    # # test_two(agent,model)
    utility_theory(agent,model)
    act!(agent, model, agent.action)    # agent performs action in agent.action
    return agent
end