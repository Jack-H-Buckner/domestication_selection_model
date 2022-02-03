"""
    HatcheryPopulations

defines a sturcture for tracking popualtions with a mix of hatchery an natural born individuals
and methods for seleting brood stock from the popuatlion. 
"""
module HatcheryPopulations


include("AgeStructuredModels.jl")

using Distributions
using DSP
using Plots
using Roots
using FFTW

mutable struct population
    natural_abundance::AbstractVector{Float64}
    hatchery_abundance::AbstractVector{Float64}
    natural_trait::AbstractMatrix{Float64}
    hatchery_trait::AbstractMatrix{Float64}
    grid::AbstractVector{Float64}
    ageStructure
    
    natural_gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    domestic_gradient::AbstractVector{Float64}
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
    correction::Float64
    
    Vle_::Float64
    s_natural::Float64
    theta_natural::Float64
    s_domestic::Float64
    theta_domestic::Float64
end 



"""
     V_star(Vle, s)

Computes equilibrium variance given severgation varinace Vle/2 and
slection strength s 
"""
function V_star(Vle, s)
    #sigma_s = 1/s
    V_prime = V -> (1/V + s)^(-1)/2 + Vle/2#V -> (V*sigma_s^2/(V + sigma_s^2))/2 + Vle/2
    V0 = Vle
    V1 = V_prime(Vle)
    while abs(V0-V1) > 10^-6
        V0 = V1
        V1 = V_prime(V0)
    end 
    return (1/V1 + s)^(-1)
end 

function init(ageStructure, Vle, theta_N, s_N, theta_D, s_D, min, max, dx)

    
    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta_N, sqrt(V_star(Vle, s_N)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),ageStructure.Amax))
    trait_A1 = transpose(repeat(transpose(trait),ageStructure.Amax))
    
    
    # Vle
    de = Distributions.Normal(theta_N, sqrt(Vle/2))
    Vle_ = Vle
    grid_Vle = collect((-3*Vle):dx:(3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    gradient_N = exp.(-s_N/2*(grid .- theta_N).^2)
    correction = 1/sum(trait .* gradient_N)
    
    gradient_D = exp.(-s_D/2*(grid .- theta_D).^2)
    
    abundance = AgeStructuredModels.stable_age_structure(ageStructure)
    #return trait
    return population(abundance, zeros(ageStructure.Amax), trait_A,trait_A1,
                    grid, ageStructure, gradient_N, gradient_D, m, Vle, correction, Vle_, 
                    s_N, theta_N, s_D, theta_D)
end 


"""
    reset!(population)

resets the popualtion to the initial caonditions specificed by init
"""
function reset!(population)
    # set trait at optimum value and variance = Vle
    d = Distributions.Normal(population.theta_natural, 
                           sqrt( V_star(population.Vle_, population.s_natural)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    
  
    d = Distributions.Normal(population.theta_natural, sqrt(V_star(population.Vle_, population.s_natural)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    trait_A1 = transpose(repeat(transpose(trait),population.ageStructure.Amax))
            
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A1
    
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
     
end 



"""
    reset!(population)

resets the popualtion to the initial caonditions specificed by init
"""
function reset!(population, theta_D)
    # set trait at optimum value and variance = Vle
    d = Distributions.Normal(population.theta_natural, 
                           sqrt( V_star(population.Vle_, population.s_natural)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    
  
    d = Distributions.Normal(population.theta_natural, sqrt(V_star(population.Vle_, population.s_natural)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    trait_A1 = transpose(repeat(transpose(trait),population.ageStructure.Amax))
            
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A1
    
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
    
    # gradient
    population.domestic_gradient = exp.(-population.s_domestic/2*(population.grid .- theta_D).^2)
    population.theta_domestic = theta_D
     
end 


"""
    reset!(population)

resets the popualtion to the initial caonditions specificed by init
"""
function reset!(population, s_N, theta_N, s_D, theta_D)
    # set trait at optimum value and variance = Vle
    
    d = Distributions.Normal(theta_N, 
                           sqrt( V_star(population.Vle_, s_N)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    
  
    d = Distributions.Normal(theta_N, sqrt(V_star(population.Vle_, s_N)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    trait_A1 = transpose(repeat(transpose(trait),population.ageStructure.Amax))
            
    population.natural_trait = trait_A
    population.hatchery_trait = trait_A1
    
    population.natural_abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
    population.hatchery_abundance = zeros(population.ageStructure.Amax)
    
    # gradient
    population.domestic_gradient = exp.(-s_D/2*(population.grid .- theta_D).^2)
    population.theta_domestic = theta_D
    population.theta_domestic = s_D
    
    population.natural_gradient = exp.(-s_N/2*(population.grid .- theta_N).^2)
    population.theta_natural = theta_N
    population.theta_natural = s_N
     
end 

end # module 