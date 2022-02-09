"""
    Methods for analyzing models that have genetic but no age structure. 
"""
module TraitStructuredPopulation


using Distributions
using DSP
using Plots
using Roots
using FFTW


mutable struct population
    ### states
    #full distribution 
    N::Float64
    trait::AbstractVector{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
    
    # normal aproximation
    mu::Float64
    V::Float64
    
    ### parameters
    # demographic
    SRCurve
    
    #full distribution
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
    correction::Float64
    
    # normal aproximation
    Vle_::Float64
    s::Float64
    theta::Float64
    
end 



mutable struct immigrants
    ### states
    # full distribuiton 
    N::Float64
    trait::AbstractVector{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
    
    # normal aproximation
    mu::Float64
    V::Float64
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


function init(SRCurve, Nstar, Vle, theta, s, min, max, dx)

    
    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    
    
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
    

    #return trait
    return population(Nstar, trait, grid, theta, V_star(Vle_, s), SRCurve, gradient, m, Vle, correction, Vle_, s, theta)
end 




function init_imigrants(population, N, mean)
    grid = population.grid
    d = Distributions.Normal(mean, sqrt(population.Vle_))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    return immigrants(N,trait,grid,mean, population.Vle_)
end 





"""
    reset!(population)

resets the popualtion to the initial caonditions specificed by init
"""
function reset!(population,Nstar)
    # set trait at optimum value and variance = Vle
    d = Distributions.Normal(population.theta, 
                           sqrt( V_star(population.Vle_, population.s)))
    trait = pdf.(d, population.grid)
    trait = trait ./sum(trait)
    population.trait = trait
    population.N = Nstar
end 


"""
    reset!(population, s)

resets the popualtion to the same initial caonditions specificed by init
and updates the value for the strength of selection 
"""
function reset!(population,s,Nstar)
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

    population.abundance = Nstar
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
    immigrants.mu = mean
end 



function reproduction(population)


    dsn = population.trait 
    

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

    return population.N, dsn
    
end 


function reproduction(population, im)

    total = population.N + im.N
    dsn = (population.N*population.trait .+ im.N*im.trait) ./ total
    

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

    return total, dsn
    
end 



"""
    reproduction_N(population)

popualtion dynamics with noraml disn and constant variance aproximation
"""
function reproduction_N(population)
    return population.N, population.mu
end 

"""
    reproduction_N(population)

popualtion dynamics with noraml disn and constant variance aproximation
"""
function reproduction_N(population, im)
    total = population.N + im.N
    p = population.N / total
    N = total*[p^2,2*p*(1-p),(1-p)^2]
    mu = [population.mu, (im.mu + population.mu)/2, im.mu, ]
    return N, mu
end 




"""
selection on juviniles 
"""
function selection(dsn, N, population)
    dsn = dsn .* population.gradient
    survival = sum(dsn)
    N = N * survival
    dsn = dsn ./ survival
    return N, dsn
end 




function recruitment(N, population)
    return population.SRCurve(N)
end 

function recruitment_N(N::AbstractVector{Float64}, population)
    Nt = population.SRCurve(sum(N))
    p = N ./ sum(N)
    return  Nt .* p
end 

function recruitment_N(N::Float64, population)
    return  population.SRCurve(N)
end 

    
function selection_N(mu::AbstractVector{Float64}, V, N::AbstractVector{Float64}, population)
    Vint = V
    V_prime = (1/Vint + population.s)^(-1)
    mu_prime = (mu./Vint .+ population.s*population.theta).*V_prime
    
    
    
    p = exp.(-1/2*((mu./sqrt(V)).^2 .+ population.s*population.theta^2 .-(mu_prime./sqrt(V_prime)).^2))
    p .*= sqrt(V_prime/V)
    
    return p.*N, mu_prime
end 

    
function selection_N(mu::Float64, V, N::Float64, population)
    Vint = V
    V_prime = (1/Vint + population.s)^(-1)
    mu_prime = (mu/Vint .+ population.s*population.theta)*V_prime
    
    
    
    p = exp(-1/2*((mu/sqrt(V))^2 + population.s*population.theta^2 -(mu_prime/sqrt(V_prime))^2))
    p *= sqrt(V_prime/V)
    

    return p*N, mu_prime
end 
    

function time_step_DSI(population)
    N,dsn = reproduction(population)
    N = recruitment(population.N, population)
    N, dsn = selection(dsn, N, population)
    population.N = N
    population.trait = dsn
end 


function time_step_DSI(population, im)
    N,dsn = reproduction(population, im)
    N = recruitment(population.N, population)
    N, dsn = selection(dsn, N, population)
    population.N = N
    population.trait = dsn
end 



function time_step_DSI_N(population)
    N,mu = reproduction_N(population)
    N = recruitment(population.N, population)
    N, mu = selection_N(mu,population.V, N, population)
    population.N = N
    population.mu = mu
end 


function time_step_DSI_N(population, im)
    N, mu = reproduction_N(population, im)
    N = recruitment_N(N, population)
    N, mu = selection_N(mu, population.V, N, population)

    population.N = sum(N)
    p = N./sum(N)
    population.mu = sum(p.*mu)
end 



function trait_moments(population)
    mu = sum(population.trait .* population.grid)
    sigma = sum(population.trait .* (population.grid .- mu).^2)
    return mu, sqrt(sigma)
end 

end # module 