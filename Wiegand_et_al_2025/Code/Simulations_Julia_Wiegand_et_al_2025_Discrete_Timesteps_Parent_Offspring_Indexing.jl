using Random, Distributions, DataFrames, Plots

# Parameters
M = 140  # Size of the torus (meters)
S = 15   # Number of species
N0 = 35  # Initial number of individuals per species
r = 8  # radius for mortality calculation
death_rate_constant = 0.05  # Constant for death rate calculation
C_effect = 0.04
H_effect = 0.04
P_L = .05
I = 0.01  # Probability of immigration event per timestep
tau = 10  # Number of timesteps to save data
Disp_k = 1.5

D_Change = 0.1  # Probability of dispersal center change per timestep

# Species-specific parameters
birth_rates = fill(0.05, S)  # Random birth rates for each species
CON = rand(S)  # Random CON values for each species
HET = rand(S)  # Random HET values for each species
K = 10  # Number of dispersal centers per species

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

# Function to calculate and plot mean CON_count and HET_count for a focal species
function plot_mean_counts(spatial_locations, focal_species)
    mean_counts = combine(groupby(spatial_locations[spatial_locations.species .== focal_species, :], :timestep), 
                          :CON_count => mean => :mean_CON_count, 
                          :HET_count => mean => :mean_HET_count)
    
    plot(mean_counts.timestep, mean_counts.mean_CON_count, label="Mean CON_count", xlabel="Timestep", ylabel="Count", title="Mean CON_count and HET_count for Species $focal_species")
    plot(mean_counts.timestep, mean_counts.mean_HET_count, label="Mean HET_count")
end

# Function to calculate and plot mean CON_count and HET_count for a focal species as a function of abundance
function plot_mean_counts_vs_abundance(spatial_locations, time_series, focal_species)
    mean_counts = combine(groupby(spatial_locations[spatial_locations.species .== focal_species, :], :timestep), 
                          :CON_count => mean => :mean_CON_count, 
                          :HET_count => mean => :mean_HET_count)
    
    abundance_data = time_series[time_series.species .== focal_species, :]
    
    merged_data = join(mean_counts, abundance_data, on=:timestep)
    
    plot(merged_data.abundance, merged_data.mean_CON_count, label="Mean CON_count", xlabel="Abundance", ylabel="Count", title="Mean CON_count vs Abundance for Species $focal_species")
    plot(merged_data.abundance, merged_data.mean_HET_count, label="Mean HET_count", xlabel="Abundance", ylabel="Count", title="Mean HET_count vs Abundance for Species $focal_species")
end



# Function to calculate and plot mean CON_count and HET_count for a focal species as a function of abundance
function plot_crowding_index_vs_abundance(spatial_locations, time_series, focal_species)
    mean_counts = combine(groupby(spatial_locations[spatial_locations.species .== focal_species, :], :timestep), 
                          :CON_count => mean => :mean_CON_count, 
                          :HET_count => mean => :mean_HET_count)
    
    abundance_data = time_series[time_series.species .== focal_species, :]
    
    merged_data = innerjoin(mean_counts, abundance_data, on=:timestep)
    
    p1 = scatter(merged_data.abundance, merged_data.mean_CON_count, label="Mean CON_count", xlabel="Abundance", ylabel="Count", title="Mean CON_count vs Abundance for Species $focal_species")
    p2 = scatter(merged_data.abundance, merged_data.mean_HET_count, label="Mean HET_count", xlabel="Abundance", ylabel="Count", title="Mean HET_count vs Abundance for Species $focal_species")
    
    plot(p1, p2, layout=(1, 2))
end

# Function to calculate and plot mean CON_count and HET_count for a focal species as a function of abundance
function plot_aggregation_metric_vs_abundance(spatial_locations, time_series, focal_species,r,M)
    mean_counts = combine(groupby(spatial_locations[spatial_locations.species .== focal_species, :], :timestep), 
                          :CON_count => mean => :mean_CON_count, 
                          :HET_count => mean => :mean_HET_count)
    
    abundance_data = time_series[time_series.species .== focal_species, :]
    
    merged_data = innerjoin(mean_counts, abundance_data, on=:timestep)
    abundance_data_others = combine(groupby(time_series[time_series.species .!= focal_species, :], :timestep), :abundance => sum => :total_abundance_others)
    abundance_vector_others = abundance_data_others.total_abundance_others

    c = 2 * π * r / M^2
    merged_data.Scaled_mean_CON_count = merged_data.mean_CON_count./(c*merged_data.abundance)
    merged_data.Scaled_mean_HET_count = merged_data.mean_HET_count./(c*abundance_vector_others)

    
    p1 = scatter(merged_data.abundance, merged_data.Scaled_mean_CON_count, label="Mean CON_count", xlabel="Abundance", ylabel="Count", title="Mean CON_count vs Abundance for Species $focal_species")
    p2 = scatter(merged_data.abundance, merged_data.Scaled_mean_HET_count, label="Mean HET_count", xlabel="Abundance", ylabel="Count", title="Mean HET_count vs Abundance for Species $focal_species")
    
    plot(p1, p2, layout=(1, 2))
end



# Initialize population
population = DataFrame(species=Int[], x=Float64[], y=Float64[], timestep=Int[], x_parent=Float64[],y_parent=Float64[],Birth=Int[])
for species in 1:S
    for _ in 1:N0
        push!(population, (species, rand() * M, rand() * M, 0,0.0,0.0,0))
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


# Function to calculate CON_count and HET_count
function calculate_counts(population, i, r, M)
    focal = population[i, :]
    neighbors = population[[toroidal_distance(focal.x, focal.y, x, y, M) <= r for (x, y) in zip(population.x, population.y)] .& (1:size(population, 1) .!= i), :]
    distances = [toroidal_distance(focal.x, focal.y, x, y, M) for (x, y) in zip(neighbors.x, neighbors.y)]
    distance_CON = distances[neighbors.species .== focal.species]
    distance_HET = distances[neighbors.species .!= focal.species]
    CON_count = isempty(distance_CON) ? 0.0 : sum(1 ./ distance_CON)
    HET_count = isempty(distance_HET) ? 0.0 : sum(1 ./ distance_HET)
    return CON_count, HET_count
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
    spatial_locations = DataFrame(species=Int[], x=Float64[], y=Float64[], x_parent=Float64[],y_parent=Float64[], Birth=Int[],timestep=Int[],CON_count=Float64[], HET_count=Float64[])

    for timestep in 1:num_timesteps
        new_population = DataFrame(species=Int[], x=Float64[], y=Float64[], timestep=Int[], x_parent=Float64[],y_parent=Float64[],Birth=Int[])
        
        # Reproduction
        for i in 1:size(population, 1)
            parent = population[i, :]
            num_offspring = rand(Poisson(birth_rates[parent.species]))
            for _ in 1:num_offspring
                (x, y) = place_new_tree(parent, dispersal_centers, P_L, M, Disp_k)
                push!(new_population, (parent.species, x, y, timestep,parent.x,parent.y,timestep))
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
            push!(new_population, (species, x, y, timestep, M/2, M, timestep))
        end

        # Dispersal center change
        for species in 1:S
            for center_index in 1:K
                if rand() < D_Change
                    dispersal_centers[species, center_index] = rand(2) * M
                end
            end
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
            #current_population = copy(population)
            #current_population.timestep .= timestep
            #append!(spatial_locations, current_population)

            for i in 1:size(population, 1)
                CON_count, HET_count = calculate_counts(population, i, r, M)
                push!(spatial_locations, (population[i, :].species, population[i, :].x, population[i, :].y,population[i, :].x_parent, population[i, :].y_parent,population[i, :].Birth, timestep, CON_count, HET_count))
            end

        end
    end

    return time_series, spatial_locations, population
end

# Run the simulation
num_timesteps = 1000  # Number of timesteps to simulate
tau = 5
Time_Series, Spatial_Locations = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)

# Call the function to plot the time series
plot_time_series(Time_Series)

# Specify the timestep you want to plot
specified_timestep = num_timesteps
# Call the function to plot the spatial locations
plot_spatial_locations(Spatial_Locations, specified_timestep)

plot_mean_counts(Spatial_Locations, 1)



# Call the updated function to plot the scaled count
focal_species = 15
plot_crowding_index_vs_abundance(Spatial_Locations, Time_Series, focal_species)

plot_aggregation_metric_vs_abundance(Spatial_Locations, Time_Series, focal_species, r, M)




plot_mean_counts_vs_abundance(Spatial_Locations,Time_Series, 1)

# Display results
println("Time Series Data:")
println(Time_Series)
println(first(Time_Series, 20))

println("Spatial Locations Data:")
println(Spatial_Locations)
println(first(Spatial_Locations, 20))



num_timesteps = 2500  # Number of timesteps to simulate
tau = 5

Random.seed!(4394)

P_L = .95
Time_Series_95, Spatial_Locations_95, Population_95 = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)

P_L = .5
Time_Series_5, Spatial_Locations_5, Population_5 = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)

P_L = .05
Time_Series_05, Spatial_Locations_05, Population_05 = run_simulation(population, dispersal_centers, birth_rates, CON, HET, C_effect, H_effect, P_L, I, r, Disp_k, D_Change, M, tau, num_timesteps)


plot_time_series(Time_Series_95)
plot_time_series(Time_Series_8)
plot_time_series(Time_Series_65)
plot_time_series(Time_Series_5)
plot_time_series(Time_Series_35)
plot_time_series(Time_Series_2)
plot_time_series(Time_Series_05)


# Call the updated function to plot the scaled count
focal_species = 15
plot_crowding_index_vs_abundance(Spatial_Locations_95, Time_Series_95, focal_species)


focal_species = 2
plot_crowding_index_vs_abundance(Spatial_Locations_05, Time_Series_05, focal_species)


plot_aggregation_metric_vs_abundance(Spatial_Locations_95, Time_Series_95, focal_species, r, M)
plot_aggregation_metric_vs_abundance(Spatial_Locations_05, Time_Series_05, focal_species, r, M)


plot_spatial_locations(Spatial_Locations_95, num_timesteps)
plot_spatial_locations(Spatial_Locations_5, num_timesteps)
plot_spatial_locations(Spatial_Locations_05, num_timesteps)

