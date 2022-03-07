module AgentBasedSimulations

include("AgeStructuredModels.jl")
include("StockRecruitCurves.jl")
using Distributions
using Roots
using StatsBase
mutable struct individual # T is an even integer represetning the ploidy of the species of interest 
    Age::Int64 # assuming age is suffincet to describe an individuals demographic state
    genes::AbstractMatrix{Float64}
    G::Float64
end 



mutable struct population
    # popualtion stat
    individuals::AbstractVector{}
    
    # demographic parameters
    ageStructureParams
    correction 
    
    # selection parameters
    s::Float64
    theta::Float64
    gradient::Function
    
    # genetic archetecture
    effect::Float64
    Ploidy::Int64 # number of chromospmes
    Nloci::Int64 # number of QTL
end 


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

function correction(Vle, s, theta, gradient)
    
    # set trait at optimum value and variance = Vle
    grid = collect(-30:0.01:30)
    d = Distributions.Normal(theta, sqrt(V_star_prime(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./ sum(trait)
    correction = 1/sum(trait .* gradient.(grid))
    
    return correction
end 



EVg(Nloci,E) = Nloci/3*E^2
#EVg(Nloci,E) = 3*Nloci*E^2
"""
age_dsn - age distribution
p - frequency at each locus 
"""
function init_individual(Amax,age_dsn, p,E,Ploidy,Nloci)
    Age = sample(1:Amax,StatsBase.Weights(age_dsn./sum(age_dsn)),1)[1]
    genes = []
    for i in 1:Nloci
        loci = []
        for j in 1:Ploidy
            if rand(1)[1] < 0.5#p[i]
                push!(loci,E/2.0)
            else
                push!(loci, -E/2.0)
            end 
        end
        push!(genes,loci)
    end 
    genes = transpose(reduce(hcat,genes))
    G = sum(genes)
  
    return individual(Age, genes, G)
end 


function init_population(ageStructureParams,s,theta,Ploidy,Nloci)
    # calcualte effect size at each locus
    E = Roots.find_zero(E -> EVg(Nloci,E) - 1.0,[0.0,10.0])
    p = rand(Nloci)
    ## compute correction used in deterministic model 
    gradient = x -> exp(-s/2*(x-theta)^2)
    correct = correction(1, s, theta, gradient)
    # get age structure 
    age_dsn = AgeStructuredModels.stable_age_structure(ageStructureParams)    
    
    Neq = floor(Int, sum(age_dsn))  
    individuals = []
    for i in 1:Neq
        push!(individuals, init_individual(ageStructureParams.Amax,age_dsn, p,E,Ploidy,Nloci))
    end 
        
    return population(individuals, ageStructureParams, correct, s,theta,gradient,E,Ploidy,Nloci)
end 

"""
    sample_genome(population)

samples genetic state of a new individual by sampling parents and then 
sampling from their allels, also allow a very small probability of mutation 

population - population object
F - precalcualted fecundities of ecah individuals 
"""
function init_juvinile(population,F)
    # sample parents
    parents = sample(population.individuals, Weights(F./sum(F)), 2, replace=true)

    # sample Haplotype 1   
    H1 = mapslices(x -> rand(x,floor(Int,population.Ploidy/2)), parents[1].genes, dims = 2)
    H2 = mapslices(x -> rand(x,floor(Int,population.Ploidy/2)), parents[2].genes, dims = 2)
    genes = hcat(H1,H2)

    # sample trait
    G = sum(genes)
    
    
    return individual(1, genes, G)
end 

function SSB(population)
    SSB = 0.0

    for ind in population.individuals
         SSB += population.ageStructureParams.Fecundity[ind.Age]
    end
    return SSB
end 
    

# before - true if selection occures before denstiy dependence 
function reproduction(population, before )
    x = SSB(population)
    N = StockRecruitCurves.density_independent(x,population.ageStructureParams.SRCurve)
    
    if before 
        N *= population.correction
    end
    
    N = rand(Distributions.Poisson(N),1)[1]
    
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    juviniles = []
    for i in 1:N
        push!(juviniles, init_juvinile(population, F))
    end 
    return juviniles
end 
using StatsBase

function reproduction_hatchery(population, N_adults, N_juviniles)
    @assert length(populations.individuals) > N_adults

    
    brood_stock = rsample(populations.individuals, N_adults, replace = false)
    
    if before 
        N_juviniles *= population.correction
    end
    
    N = rand(Distributions.Poisson(N),1)[1]
    
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    juviniles = []
    for i in 1:N
        push!(juviniles, init_juvinile(population, F))
    end 
    return juviniles
end 

function selection(population, juviniles)
    survival = zeros(length(juviniles))
    i=0
    for ind in juviniles
        i+=1
        if population.gradient(ind.G) > rand(1)[1]
            survival[i] = 1
        else
            survival[i] = 0
        end
    end 
    return juviniles[survival.==1]
end 
    
    
function desnity_dependence(population, juviniles, before)
    N = length(juviniles)
    if before
        p = StockRecruitCurves.density_dependent(N,population.ageStructureParams.SRCurve)
        N = rand(Distributions.Binomial(N,p), 1)[1]
    else
        p = StockRecruitCurves.density_dependent(N,population.ageStructureParams.SRCurve)
    
        p*=population.correction
        N = rand(Distributions.Binomial(N,p), 1)[1]
    end 
    return juviniles[1:N]
end


function directional_immigration!(p,juviniles_in, juviniles_out)
    N = rand(Distributions.Binomial(length(juviniles_out),p))

    if N >0
        immigrants = juviniles_out[1:N]
        juviniles_out = juviniles_out[N:end]
        juviniles_in = vcat(juviniles_in, immigrants)
    end
    return juviniles_in, juviniles_out
end 
    
function immigration!(p,juviniles1, juviniles2)
    N1 = rand(Distributions.Binomial(length(juviniles1),p))
    N2 = rand(Distributions.Binomial(length(juviniles2),p))
    if N1 > 0
        immigrants1 = juviniles1[1:N1]
        juviniles1 = juviniles1[N1:end]
        juviniles2 =vcat(juviniles2,immigrants1)
    end 
    if N2 > 0
        immigrants2 = juviniles1[1:N2]
        juviniles2 = juviniles2[N2:end]
        juviniles1 = vcat(juviniles1,immigrants2)
    end 
    return juviniles1, juviniles2
end
    
    
function update_population!(population, juviniles)
    # survival + aging
    survival = zeros(length(population.individuals))
    psurvival = population.ageStructureParams.Survival
    psurvival[end] = 0
    i = 0
    for ind in population.individuals
        i+=1
        
        if ind.Age == 150
            survival[i] = 0.0
        elseif ind.Age > 150
            survival[i] = 0.0
        else
            survival[i] = 1 * (psurvival[ind.Age] > rand(1)[1])
        end
        ind.Age+=1
    end 
    
    population.individuals = population.individuals[survival .==1]
        
    # add juviniles 
    for ind in juviniles
        push!(population.individuals, ind)
    end 
    #population.individuals = population.individuals[broadcast(x -> x.Age, population.individuals) .< 151]
    
end 
    
    
### single population dynamics
    
function time_step!(population, before)
    if before
        J = reproduction(population, before )
        J = selection(population, J)
        J = desnity_dependence(population, J, before)
        update_population!(population, J)
    else
        J = reproduction(population, before )
        J = desnity_dependence(population, J, before)
        J = selection(population, J)
        update_population!(population, J)
    end
end 


### single population dynamics
    
function time_step!(population1, population2, p_sym, p_asym)

        J1 = reproduction(population1, false )
        J1 = desnity_dependence(population1, J1, false)
        J1 = selection(population1, J1)
        
        
        J2 = reproduction(population2, false )
        J2 = desnity_dependence(population2, J2, false)
        J2 = selection(population2, J2)
        
        
        J1, J2 = immigration!(p_sym,J1, J2)
        J1, J2 = directional_immigration!(p_asym,J1, J2)
        
#         maxAge = broadcast(x->x.Age, population1.individuals)[argmax(broadcast(x->x.Age, population1.individuals))]
#     maxAge2 = broadcast(x->x.Age, population1.individuals)[argmax(broadcast(x->x.Age, population1.individuals))]
#     print(maxAge2)
#     print(" ")
#         println(maxAge)
        update_population!(population1, J1)
        update_population!(population2, J2)
end 



## population statistics

function trait_distribution(population)
    Gls = []
    for ind in population.individuals
        push!(Gls, ind.G)
    end 
    return Gls 
end 

function trait_distribution_moments(population)
    Gls = trait_distribution(population)
    mu = sum(Gls)/length(Gls)
    sigma = sqrt(sum((Gls .- mu).^2 ./(length(Gls) - 1)))
    return mu, sigma
end 



function trait_distribution_age1(population)
    Gls = []
    for ind in population.individuals
        if ind.Age == 1
            push!(Gls, ind.G)
        end
    end 
    return Gls 
end 

function trait_distribution_moments_age1(population)
    Gls = trait_distribution_age1(population)
    mu = sum(Gls)/length(Gls)
    sigma = sqrt(sum((Gls .- mu).^2 ./(length(Gls) - 1)))
    return mu, sigma
end 




function age_distribution(population)
    Gls = []
    for ind in population.individuals
        push!(Gls, ind.Age)
    end 
    return Gls 
end 

end # module 