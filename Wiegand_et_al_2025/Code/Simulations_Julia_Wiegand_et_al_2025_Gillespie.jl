using Random, Distributions, DataFrames, Plots

# Parameters
M = 100  # Size of the torus (meters)
S = 10    # Number of species
N0 = 75  # Initial number of individuals per species
r = 10  # Radius for mortality calculation
C_effect = 0.0075
H_effect = 0.0075
P_L = .05
I = 0.01  # Immigration rate
tau = 100  # Save data every tau events
Disp_k = 5
D_Change = 0.01  # Rate at which a dispersal center changes location

# Species-specific parameters
birth_rates = fill(0.2, S)  # Random birth rates for each species
CON = rand(S)  # Random CON values for each species
HET = rand(S)  # Random HET values for each species
K = 25  # Number of dispersal centers per species

# Initialize dispersal centers
dispersal_centers = [rand(2) * M for _ in 1:S, _ in 1:K]

# Plot the time series
function plot_time_series(time_series)
    species_list = unique(time_series.species)
    plot()
    for species in species_list
        species_data = time_series[time_series.species .== species, :]
        plot!(species_data.event, species_data.abundance, label="Species $species")
    end
    xlabel!("Event Count")
    ylabel!("Abundance")
    title!("Abundance of Each Species Through Time")
end

# Plot spatial locations
function plot_spatial_locations(spatial_locations, event_count)
    event_data = spatial_locations[spatial_locations.event .== event_count, :]
    species_list = unique(event_data.species)
    plot()
    for species in species_list
        species_data = event_data[event_data.species .== species, :]
        scatter!(species_data.x, species_data.y, label="Species $species")
    end
    xlabel!("X Coordinate")
    ylabel!("Y Coordinate")
    title!("Spatial Locations of Individuals at Event $event_count")
end

# Initialize population
population = DataFrame(species=Int[], x=Float64[], y=Float64[], event=Int[], time=Float64[])
for species in 1:S
    for _ in 1:N0
        push!(population, (species, rand() * M, rand() * M, 0, 0.0))
    end
end

# Function to calculate toroidal distance
function toroidal_distance(x1, y1, x2, y2, M)
    dx = abs(x1 - x2)
    dy = abs(y1 - y2)
    dx = min(dx, M - dx)
    dy = min(dy, M - dy)
    return sqrt(dx^2 + dy^2)
end

# Function to calculate mortality rate
function mortality_rate(population, i, r, C_effect, H_effect, CON, HET, M)
    focal = population[i, :]
    neighbors = population[[toroidal_distance(focal.x, focal.y, x, y, M) <= r for (x, y) in zip(population.x, population.y)], :]
    CON_count = sum(neighbors.species .== focal.species)
    HET_count = size(neighbors, 1) - CON_count
    return CON_count * C_effect + HET_count * H_effect
end

# Function to place a new tree
function place_new_tree(parent, dispersal_centers, P_L, M, Disp_k)
    if rand() < P_L
        # Place near parent
        x = mod(parent.x + rand(Exponential(Disp_k)), M)
        y = mod(parent.y + rand(Exponential(Disp_k)), M)
    else
        # Place near dispersal center
        center = dispersal_centers[parent.species, rand(1:K)]
        x = mod(center[1] + rand(Exponential(Disp_k)), M)
        y = mod(center[2] + rand(Exponential(Disp_k)), M)
    end
    return (x, y)
end

# Gillespie algorithm function
function run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_events)
    event_count = 0
    absolute_time = 0.0
    time_series = DataFrame(event=Int[], time=Float64[], species=Int[], abundance=Int[])
    spatial_locations = DataFrame(species=Int[], x=Float64[], y=Float64[], event=Int[], time=Float64[])

    while event_count < num_events
        # Calculate rates
        death_rates = [mortality_rate(population, i, r, C_effect, H_effect, CON, HET, M) for i in 1:size(population, 1)]
        birth_rates_individual = [birth_rates[population[i, :].species] for i in 1:size(population, 1)]
        total_rate = sum(death_rates) + sum(birth_rates_individual) + I + D_Change

        # Determine time to next event
        dt = rand(Exponential(1 / total_rate))
        absolute_time += dt
        event_count += 1

        # Determine which event occurs
        event_prob = rand() * total_rate
        if event_prob < sum(death_rates)
            # Death event
            i = findfirst(cumsum(death_rates) .>= event_prob)
            delete!(population, i)
        elseif event_prob < sum(death_rates) + sum(birth_rates_individual)
            # Birth event
            i = findfirst(cumsum(birth_rates_individual) .>= event_prob - sum(death_rates))
            parent = population[i, :]
            (x, y) = place_new_tree(parent, dispersal_centers, P_L, M, Disp_k)
            push!(population, (parent.species, x, y, event_count, absolute_time))
        elseif event_prob < sum(death_rates) + sum(birth_rates_individual) + I
            # Immigration event
            species = rand(1:S)
            x = rand() * M
            y = rand() * M
            push!(population, (species, x, y, event_count, absolute_time))
        else
            # Dispersal center change event
            species = rand(1:S)
            center_index = rand(1:K)
            dispersal_centers[species, center_index] = rand(2) * M
        end

        # Save data every tau events
        if event_count % tau == 0
      #      println("Event $event_count: Saving data...")
            # Save abundance data
            abundance = combine(groupby(population, :species), nrow => :count)
            for row in eachrow(abundance)
                push!(time_series, (event_count, absolute_time, row.species, row.count))
            end
            # Save spatial data with the current event_count and absolute_time
            current_population = copy(population)
            current_population.event .= event_count
            current_population.time .= absolute_time
            append!(spatial_locations, current_population)
        end
    end

    return time_series, spatial_locations
end

# Run the simulation
num_events = 500  # Number of events to simulate
tau = 25
Time_Series, Spatial_Locations = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_events)

# Display results
println("Time Series Data:")
println(Time_Series)
println(first(Time_Series, 20))

println("Spatial Locations Data:")
println(Spatial_Locations)
println(first(Spatial_Locations, 20))

# Call the function to plot the time series
plot_time_series(Time_Series)

# Specify the event count you want to plot
specified_event_count = 5000
# Call the function to plot the spatial locations
plot_spatial_locations(Spatial_Locations, specified_event_count)