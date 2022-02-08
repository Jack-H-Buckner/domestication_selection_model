"""
ths file provides funciton that define an age structured popualtion model 
that tracks he full age dependent distribution of a trait in the population. 

The code is built around a centeral mutable struct, 'population'. Popualtion stores 
a vector with abundcaes of each age class, and an array that stores the density of the 
trait distribution for each cohort at a set of nodes. In additon, to these states the object
stores paramters that describe the fecudndity and survival of each age class, a stock - recruit
relationship, selection gradient and trait variance. 

Around this basic structure, several methods are defined to update the state of the popuatlion
as as differnt processes occur, such as immigraiton, reproduction, aging and recruitment. 
"""
module AgeTraitStructuredModels

include("AgeStructuredModels.jl")

using Distributions
using DSP
using Plots
using Roots
using FFTW



mutable struct population
    # states
    abundance::AbstractVector{Float64}
    trait::AbstractMatrix{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
    
    # parameters
    ageStructure
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
    correction::Float64
    
    Vle_::Float64
    s::Float64
    theta::Float64
end 


mutable struct immigrants
    # states
    N::Float64
    trait::AbstractVector{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
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

"""
     V_star_prime(Vle, s)

Computes equilibrium variance given severgation varinace Vle/2 and
slection strength s. measurement after reproduciton but before selection
"""
function V_star_prime(Vle, s)
    #sigma_s = 1/s
    V_prime = V -> (1/V + s)^(-1)/2 + Vle/2#V -> (V*sigma_s^2/(V + sigma_s^2))/2 + Vle/2
    V0 = Vle
    V1 = V_prime(Vle)
    while abs(V0-V1) > 10^-6
        V0 = V1
        V1 = V_prime(V0)
    end 
    return V1
end 

function init(ageStructure, Vle, theta, s, min, max, dx)

    
    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),ageStructure.Amax))
    
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    Vle_ = Vle
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    d = Distributions.Normal(theta, sqrt(V_star_prime(Vle_, s)))
    trait = pdf.(d, grid)
    trait = trait ./ sum(trait)
    gradient = exp.(-s/2*(grid .- theta).^2)
    correction = 1/sum(trait .* gradient)
    



    abundance = AgeStructuredModels.stable_age_structure(ageStructure)
    #return trait
    return population(abundance, trait_A, grid, ageStructure, gradient, m, Vle, correction, Vle_, s, theta)
end 




function init_imigrants(population, N, mean)
    grid = population.grid
    d = Distributions.Normal(mean, sqrt(population.Vle_))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    return immigrants(N,trait,grid)
end 


"""
    reset!(population)

resets the popualtion to the initial caonditions specificed by init
"""
function reset!(population)
    # set trait at optimum value and variance = Vle
    d = Distributions.Normal(population.theta, 
                           sqrt( V_star(population.Vle_, population.s)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    population.trait = trait_A
    populations.abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
end 


"""
    reset!(population, s)

resets the popualtion to the same initial caonditions specificed by init
and updates the value for the strength of selection 
"""
function reset!(population,s)
    # set trait at optimum value and variance = Vle
    d = Distributions.Normal(population.theta, 
                            sqrt(V_star(population.Vle_,s)))
    #println(sqrt(V_star(population.Vle_,s)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),population.ageStructure.Amax))
    population.trait = trait_A
    
    # gradient
    d = Distributions.Normal(population.theta, 
                            sqrt(V_star_prime(population.Vle_,s)))
    trait = pdf.(d, population.grid)
    trait = trait ./ sum(trait)
    trait = trait ./sum(trait)
    
    population.gradient = exp.(-s/2*(population.grid .- population.theta).^2)
    population.correction =  1/sum(trait .* population.gradient)

    population.abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
end 

"""
    reset_immigrants!(im, N, mean)

updates immigants object with new mena vlaue and m
"""
function reset_immigrants!(immigrants,population, N, mean)
    d = Distributions.Normal(mean, sqrt(population.Vle_))
    immigrants.trait = Distributions.pdf.(d, population.grid)
    immigrants.trait= immigrants.trait ./sum(immigrants.trait)
    immigrants.N = N
end 
"""
returns trait distribution for new age class and 
the total spawning stock fecundity.
"""
function reproduction_fft(population)

    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    dsn = population.trait * (population.abundance .* 
                population.ageStructure.Fecundity) ./ f_total
    

    #Plots.plot!(population.grid,dsn )
    # convolution - random mating 
    N = length(population.grid)-1
    dsn_zs = vcat(dsn,zeros(length(dsn)))
    dsn = real(ifft(fft(dsn_zs).*fft(dsn_zs)))[1:2:(2*N+1)]#DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    dsn = dsn./sum(dsn)
    

    #Plots.plot!(population.grid,dsn )

    # convolution inperfect inheritance 
    m = convert(Int64,floor(population.m/2+1))
    dsn = DSP.conv(dsn, population.Vle)
 
    dsn= dsn[m:(N+m)]
    dsn = dsn ./ sum(dsn)

 
    return dsn, f_total
end 

function reproduction(population)

    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    dsn = population.trait * (population.abundance .* 
                                population.ageStructure.Fecundity) ./ f_total

    #Plots.plot!(population.grid,dsn )
    # convolution - random mating 
    N = length(population.grid)-1
    dsn = DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    dsn = dsn./sum(dsn)
    

    #Plots.plot!(population.grid,dsn )

    # convolution inperfect inheritance 
    m = convert(Int64,floor(population.m/2+1))
    dsn = DSP.conv(dsn, population.Vle)
 
    dsn= dsn[m:(N+m)]
    dsn = dsn ./ sum(dsn)

 
    return dsn, f_total
end 
"""
selection on juviniles 
"""
function selection(dsn, N, population)
    dsn = dsn .* population.gradient
    survival = sum(dsn)
    N = N * survival
    dsn = dsn ./ survival
    return dsn, N
end 



"""
Beverton Holt recruitment 
assumes that selection occurs after recruitment 
"""
function recruitment(f_total, population )
    R = population.correction * population.ageStructure.SRCurve(f_total) #
    return R 
end 

"""
Beverton Holt recruitment 

f_total - number of eggs
populations - populaiton objct for parameters
before - boolian does selection occur before or after selection 
"""
function recruitment(f_total, population, before)
    if before 
        f_total *= population.correction 
        R = population.ageStructure.SRCurve(f_total)
    else
        R = population.correction * population.ageStructure.SRCurve(f_total)
    end 
    return R 
end 



"""
    immigration(dsn,N, immigrants)

Immigration of juviniles. This function is designed for populations 
"""
function immigration(dsn,N, immigrants)
    N_total = N + immigrants.N
    p = N/N_total
    dsn = p.*dsn .+ (1-p).* immigrants.trait
    return dsn, N_total
end 




"""
returns trait distribution for new age class and 
the total spawning stock fecundity.
"""
function selection_and_reproduction(population)
    
    S = transpose(repeat(transpose(population.gradient),population.A_max))
    
    dsn = population.trait .* S * (population.abundance .*  population.fecundity) 
    f_total = sum(dsn)
    dsn = dsn./f_total
    
    
    # convolution - random mating 
    N = length(population.grid)-1
    dsn = DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    dsn = dsn./sum(dsn)
    
    # convolution inperfect inheritance 
    m = convert(Int64,floor(population.m/2+1))
    dsn = DSP.conv(dsn, population.Vle)[m:(N+m)]
    dsn = dsn ./ sum(dsn)
    
    return dsn, f_total
end 




 

"""
Updates the age structure of the popuatlion and adds recruits
"""
function ageing!(population, R, dsn_R)
    population.abundance = population.abundance .* population.ageStructure.Survival
    plus_group = sum(population.abundance[end-1:end])
    new_A = zeros(population.ageStructure.Amax)
    new_A[1] = R 
    #new_A[end] = plus_group
    

    
    new_A[2:end] = population.abundance[1:end-1]
    
    N = length(population.grid)
    new_dsn = zeros(N,population.ageStructure.Amax)
    new_dsn[:,1] = dsn_R
    #plus_group_dsn = (population.abundance[end-1]*population.trait[:,end-1] .+ population.abundance[end]*population.trait[:,end]) ./plus_group
    #new_dsn[:,end] = plus_group_dsn
    
    new_dsn[:,2:end] = population.trait[:,1:end-1]
    
    population.abundance = new_A
    population.trait = new_dsn
    
    #return population
end 

"""
    simulates dynamics and updates the correction factor once trait dsn has equilibrated 
"""
function update_correction!(population)
    for i in 1:1000
        dsn, R = reproduction(population)
        R = recruitment(R, population )
        dsn, R = selection(dsn, R, population)
        ageing!(population, R, dsn)
    end 
    dsn, R = reproduction(population)
    dsn, R = selection(dsn, 1, population)
    population.correction = 1/R
    population.abundance = AgeStructuredModels.stable_age_structure(population.ageStructure)
end 

###################################
### time step/ update functions ###
###################################

"""
    time_step_DSI!(populations; immigrants)
updates the popualtion with density dependnece before selection. If immigrants are included
immigration occurs last
"""
function time_step_DSI!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_DSI!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


"""
    time_step_DIS!(populations; immigrants)
updates the popualtion with density dependnece before selection. If immigrants are included
immigration occurs before selection.
"""
function time_step_DIS!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_DIS!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 



"""
    time_step_IDS!(populations; immigrants)
updates the popualtion with density dependnece before selection. If immigrants are included
immigration occurs before density dependence.
"""
function time_step_IDS!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_IDS!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    R = AgeTraitStructuredModels.recruitment(R, population )
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 



"""
    time_step_ISD!(populations; immigrants)
updates the popualtion with selection before density dependnece. If immigrants are included
immigration occurs before selection.
"""
function time_step_ISD!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_ISD!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


"""
    time_step_ISD!(populations; immigrants)
updates the popualtion with selection before density dependnece. If immigrants are included
immigration occurs before selection.
"""
function time_step_SDI!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_SDI!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


"""
    time_step_ISD!(populations; immigrants)
updates the popualtion with selection before density dependnece. If immigrants are included
immigration occurs before selection.
"""
function time_step_SID!(population)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 


function time_step_SID!(population, immigrants)
    dsn, R = AgeTraitStructuredModels.reproduction(population)
    dsn, R = AgeTraitStructuredModels.selection(dsn, R, population)
    dsn, R = AgeTraitStructuredModels.immigration(dsn, R, immigrants)
    R = AgeTraitStructuredModels.recruitment(R, population, true )
    AgeTraitStructuredModels.ageing!(population, R, dsn)
end 





## analysis
"""
This funcction computes fittness in terms of the average 
probability of survival of juviniles produved by a popualtion
assuming selection acts on the survival of juviniles to recruitment. 
"""
function fittness(population)
    dsn, f = reproduction(population)
    W = sum(dsn .* population.gradient)
    return W
end 


# 
## analysis
"""
computes the number of individuals of age 1 in the population
"""
function recruitment(population)
    return population.abundance[1]
end 


"""
computes the number of individuals of age 1 in the population
"""
function spawning_stock(population)
    max_F = argmax(population.ageStructure.Fecundity)
    return sum(population.abundance .* population.ageStructure.Fecundity) ./ population.ageStructure.Fecundity[max_F]
end 


"""
    trait_moments_spawning(population)

returns the first two moments of the spawning stocks trait distirbution
"""
function trait_moments_spawning(population)
    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    dsn = population.trait * (population.abundance .* 
                population.ageStructure.Fecundity) ./ f_total
    mu = sum(dsn .* population.grid)
    sigma = sum(dsn .* (population.grid .- mu).^2)
    return mu, sqrt(sigma)
end 


"""
    trait_moments_spawning(population)

returns the first two moments of the spawning stocks trait distirbution
"""
function trait_moments_born(population)
    dsn, R = reproduction(population)
    mu = sum(dsn .* population.grid)
    sigma = sum(dsn .* (population.grid .- mu).^2)
    return mu, sqrt(sigma)
end 

"""
    trait_moments_recruits(population)

returns the first two moments of the age one trait distirbution
"""
function trait_moments_recruits(population)
    dsn = population.trait[:,1] ./ sum(population.trait[:,1] )
    mu = sum(dsn .* population.grid)
    sigma = sum(dsn .* (population.grid .- mu).^2)
    return mu, sqrt(sigma)
end 

end # module 




#function init_population(A_max, survival, fecundity, r , K, theta, s, min, max, dx, Vle)
    
#     # set age distribution - equilibrium 
#     LEP = age_structure_model.LEP(1.0, survival, fecundity, A_max)
    
#     a = r/LEP
#     b = a/K
#     f(x) = a*x/(1+b*x)
#     g(x )= f(x) - x/LEP
#     print(LEP)
#     x_eq = Roots.find_zero(g, (10^-6, 10^6*K*LEP))
#     y_eq = f(x_eq)
    
#     abundance = zeros(A_max)
#     N = y_eq
#     for i in 1:A_max
#         abundance[i] = N
#         N *= survival[i]
#     end 

#     # set trait at optimum value and variance = Vle
#     grid = collect(min:dx:max)
#     d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
#     trait = pdf.(d, grid)
#     trait = trait ./sum(trait)
#     trait_A = transpose(repeat(transpose(trait),A_max))
    
    
#     # Vle
#     de = Distributions.Normal(theta, sqrt(Vle/2))
#     grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
#     Vle = pdf.(de,grid_Vle)
#     m = length(Vle)
    
#     # gradient
#     gradient = exp.(-s/2*(grid .- theta).^2)
#     correction = 1/sum(trait .* gradient)
    
#     pop = population(abundance, trait_A, grid, fecundity, survival,A_max, a, b, correction, gradient,m,Vle)
    
#     return pop
# end 



# function init_population_2(A_max, survival, fecundity, R_star, alpha , beta, theta, s, min, max, dx, Vle, correction)
    
#     # set age distribution - equilibrium 
#     LEP = age_structure_model.LEP(1.0, survival, fecundity, A_max)
    
    
#     abundance = zeros(A_max)
#     N = R_star
#     for i in 1:A_max
#         abundance[i] = N 
#         N *= survival[i]
#     end 
#     # set trait at optimum value and variance = Vle 
#     grid = collect(min:dx:max)
#     d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
#     trait = pdf.(d, grid)
#     trait = trait ./sum(trait)
#     trait_A = transpose(repeat(transpose(trait),A_max))
    
    
#     # Vle
#     de = Distributions.Normal(theta, sqrt(Vle/2))
#     grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
#     Vle = pdf.(de,grid_Vle)
#     m = length(Vle)
    
#     # gradient
#     if correction
#         gradient = exp.(-s/2*(grid .- theta).^2)
#         correction = 1/sum(trait .* gradient)
#     else
#         correction = 1
#     end 
    
    
    
#     pop = population(abundance, trait_A, grid, fecundity, survival,A_max, alpha, beta, correction, gradient,m,Vle)

#     return pop
# end 
