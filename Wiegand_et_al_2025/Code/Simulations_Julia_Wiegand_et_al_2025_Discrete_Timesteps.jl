using Random, Distributions, DataFrames, Plots

# Parameters
M = 100  # Size of the torus (meters)
S = 10   # Number of species
N0 = 15  # Initial number of individuals per species
r = 0.1  # Average reproduction rate per individual per timestep
death_rate_constant = 0.05  # Constant for death rate calculation
C_effect = 0.05
H_effect = 0.05
P_L = 1
I = 0.01  # Probability of immigration event per timestep
tau = 10  # Number of timesteps to save data
Disp_k = 5
D_Change = 0.01  # Probability of dispersal center change per timestep

# Species-specific parameters
birth_rates = fill(0.3, S)  # Random birth rates for each species
CON = rand(S)  # Random CON values for each species
HET = rand(S)  # Random HET values for each species
K = 5  # Number of dispersal centers per species

# Initialize dispersal centers
dispersal_centers = [rand(2) * M for _ in 1:S, _ in 1:K]

# Plot the time series
function plot_time_series(time_series)
    species_list = unique(time_series.species)
    plot()
    for species in species_list
        species_data = time_series[time_series.species .== species, :]
        plot!(species_data.timestep, species_data.abundance, label="Species $species")
    end
    xlabel!("Timestep")
    ylabel!("Abundance")
    title!("Abundance of Each Species Through Time")
end

# Plot spatial locations
function plot_spatial_locations(spatial_locations, timestep)
    timestep_data = spatial_locations[spatial_locations.timestep .== timestep, :]
    species_list = unique(timestep_data.species)
    plot()
    for species in species_list
        species_data = timestep_data[timestep_data.species .== species, :]
        scatter!(species_data.x, species_data.y, label="Species $species")
    end
    xlabel!("X Coordinate")
    ylabel!("Y Coordinate")
    title!("Spatial Locations of Individuals at Timestep $timestep")
end

# Initialize population
population = DataFrame(species=Int[], x=Float64[], y=Float64[], timestep=Int[])
for species in 1:S
    for _ in 1:N0
        push!(population, (species, rand() * M, rand() * M, 0))
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
    neighbors = population[[toroidal_distance(focal.x, focal.y, x, y, M) <= r for (x, y) in zip(population.x, population.y)] .& (1:size(population, 1) .!= i), :]
    distances = [toroidal_distance(focal.x, focal.y, x, y, M) for (x, y) in zip(neighbors.x, neighbors.y)]
    distance_CON = distances[neighbors.species .== focal.species]
    distance_HET = distances[neighbors.species .!= focal.species]
    CON_count = isempty(distance_CON) ? 0.0 : sum(1 ./ distance_CON)
    HET_count = isempty(distance_HET) ? 0.0 : sum(1 ./ distance_HET)
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

# Simulation function
function run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)
    time_series = DataFrame(timestep=Int[], species=Int[], abundance=Int[])
    spatial_locations = DataFrame(species=Int[], x=Float64[], y=Float64[], timestep=Int[])

    for timestep in 1:num_timesteps
        new_population = DataFrame(species=Int[], x=Float64[], y=Float64[], timestep=Int[])

        # Reproduction
        for i in 1:size(population, 1)
            parent = population[i, :]
            num_offspring = rand(Poisson(r))
            for _ in 1:num_offspring
                (x, y) = place_new_tree(parent, dispersal_centers, P_L, M, Disp_k)
                push!(new_population, (parent.species, x, y, timestep))
            end
        end

        # Death
        death_rates = [mortality_rate(population, i, r, C_effect, H_effect, CON, HET, M) for i in 1:size(population, 1)]
        survivors = population[[rand() < exp(-death_rate) for death_rate in death_rates], :]
        append!(new_population, survivors)

        #print(survivors)
        #println(new_population)

        # Immigration
        if rand() < I
            species = rand(1:S)
            x = rand() * M
            y = rand() * M
            push!(new_population, (species, x, y, timestep))
        end

        # Dispersal center change
        if rand() < D_Change
            species = rand(1:S)
            center_index = rand(1:K)
            dispersal_centers[species, center_index] = rand(2) * M
        end

        population = new_population

        # Save data every tau timesteps
        if timestep % tau == 0
            println("Timestep $timestep: Saving data...")
            # Save abundance data
            abundance = combine(groupby(population, :species), nrow => :count)
            for row in eachrow(abundance)
                push!(time_series, (timestep, row.species, row.count))
            end
            # Save spatial data with the current timestep
            append!(spatial_locations, population)
        end
    end

    return time_series, spatial_locations
end

# Run the simulation
num_timesteps = 1000  # Number of timesteps to simulate
tau = 10
Time_Series, Spatial_Locations = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)

# Display results
println("Time Series Data:")
println(Time_Series)
println(first(Time_Series, 20))

println("Spatial Locations Data:")
println(Spatial_Locations)
println(first(Spatial_Locations, 20))

# Call the function to plot the time series
plot_time_series(Time_Series)

# Specify the timestep you want to plot
specified_timestep = 50
# Call the function to plot the spatial locations
plot_spatial_locations(Spatial_Locations, specified_timestep)