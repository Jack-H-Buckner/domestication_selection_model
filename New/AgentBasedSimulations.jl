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
    ID::Float64
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


function V_star(Vle, s)
    #sigma_s = 1/s
    V_prime = V -> (1/V + s)^(-1)/2 + Vle/2#V -> (V*sigma_s^2/(V + sigma_s^2))/2 + Vle/2
    V0 = Vle
    V1 = V_prime(Vle)
    while abs(V0-V1) > 10^-6
        V0 = V1
        V1 = V_prime(V0)
    end 
    return (1/V1 + s)^(-1)/2 + Vle/2
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


function WF_probabilities(p,ploidy)
    d = Distributions.Binomial(ploidy,p)
    return pdf.(d, 0:ploidy)
end 

function WF_levels(E,ploidy)
    levels = 2*E .* collect(0:ploidy)
    levels .+= -1*E*ploidy
    return levels
end 

function WF_variance(p,E,ploidy)
    mu = sum(WF_levels(E,ploidy) .* WF_probabilities(p,ploidy))
        
    return sum(WF_probabilities(p,ploidy).*(WF_levels(E,ploidy) .- mu).^2 ) 
end 

function WF_expected_variance(E,ploidy)
    p = 0:0.005:1
    
    return sum(WF_variance.(p,E,ploidy) .* 0.005)
end 

function WF_trait_variance(E, Nloci,ploidy)
    return Nloci*WF_expected_variance(E,ploidy)
end 


function solve_WF_effects(variance, Nloci,ploidy)
    return Roots.find_zero(E -> WF_trait_variance(E, Nloci,ploidy) - variance, [0.0, 2*variance])
end 

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
            if rand(1)[1] < p[i]
                push!(loci,E)
            else
                push!(loci, -E)
            end 
        end
        push!(genes,loci)
    end 
    genes = transpose(reduce(hcat,genes))
    
    G = sum(genes)
  
    return individual(Age, genes, G, rand(1)[1])
end 


function init_population(ageStructureParams,s,theta,Ploidy,Nloci)
    # calcualte effect size at each locus
    E = solve_WF_effects(1.0, Nloci,Ploidy)
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
    
    
    return individual(1, genes, G, rand(1)[1])
end 

function SSB(population)
    SSB = 0.0

    for ind in population.individuals
        if ind.Age > 150
            print("here")
        end
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
     
        
function reproduction(population, before, sigma )
    x = SSB(population)
    N = StockRecruitCurves.density_independent(x,population.ageStructureParams.SRCurve)
    d = Distributions.Normal(-1*(sigma^2)/2,sigma)
    if before 
        N *= population.correction*exp(rand(d))
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

        
"""
samples a sub set of adults to use as brood stock 
"""
function reproduction_hatchery(population, N_adults, N_juviniles, before)
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
        if (p > 0) & (p < 1)
            N = rand(Distributions.Binomial(N,p), 1)[1]           
        else
            N = 0
        end
    end 
    return juviniles[1:N]
end


function directional_immigration!(p,juviniles_in, juviniles_out)
    N = rand(Distributions.Binomial(length(juviniles_out),p))
    if N >0
        immigrants = juviniles_out[1:N]
        juviniles_out = juviniles_out[(N+1):end]
        juviniles_in = vcat(juviniles_in, immigrants)
    end   
    return juviniles_in, juviniles_out
end 
    
function immigration!(p,juviniles1, juviniles2)
    N1 = rand(Distributions.Binomial(length(juviniles1),p))
        
    N2 = rand(Distributions.Binomial(length(juviniles2),p))

                
    if N1 > 0
        immigrants1 = juviniles1[1:N1]
        juviniles1 = juviniles1[(N1+1):end]
        juviniles2 =vcat(juviniles2,immigrants1)
    end 
    if N2 > 0
        immigrants2 = juviniles2[1:N2]
        juviniles2 = juviniles2[(N2+1):end]
        juviniles1 = vcat(juviniles1,immigrants2)
    end 
#     print(" ")
#     println(length(juviniles1) + length(juviniles2))
            
    return juviniles1, juviniles2
end
    
    
function update_population!(population, juviniles)
    # survival + aging
    survival = zeros(length(population.individuals))
    psurvival = population.ageStructureParams.Survival
    psurvival[end] = 0

    for i in 1:length(population.individuals) 
            
        ind = population.individuals[i]
        if ind.Age > 149
            survival[i] = 0.0
        else
            survival[i] = 1.0 * (psurvival[ind.Age] > rand(1)[1])
        end
        ind.Age+=1
    end 
   
    
    population.individuals = population.individuals[survival .==1.0]
    
  
    
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
        J1 = desnity_dependence(population1,J1,false)
        J1 = selection(population1,J1)
        
        
        J2 = reproduction(population2,false)
        J2 = desnity_dependence(population2,J2,false)
        J2 = selection(population2,J2)
        
        
        J1, J2 = immigration!(p_sym,J1,J2)
        J1, J2 = directional_immigration!(p_asym,J1, J2)
        
        update_population!(population1, J1)
        update_population!(population2, J2)
end 


    
function time_step!(population, before, sigma)
    if before
        J = reproduction(population, before, sigma )
        J = selection(population, J)
        J = desnity_dependence(population, J, before)
        update_population!(population, J)
    else
        J = reproduction(population, before, sigma )
        J = desnity_dependence(population, J, before)
        J = selection(population, J)
        update_population!(population, J)
    end
end 


### single population dynamics
    
function time_step!(population1, population2, p_sym, p_asym, sigma)

        J1 = reproduction(population1, false, sigma )
        J1 = desnity_dependence(population1,J1,false)
        J1 = selection(population1,J1)
        
        
        J2 = reproduction(population2,false, sigma)
        J2 = desnity_dependence(population2,J2,false)
        J2 = selection(population2,J2)
        
        
        J1, J2 = immigration!(p_sym,J1,J2)
        J1, J2 = directional_immigration!(p_asym,J1, J2)
        
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
    if length(Gls) > 0 
        mu = sum(Gls)/length(Gls)
        sigma = sqrt(sum((Gls .- mu).^2 ./(length(Gls) - 1)))
        return mu, sigma
    else
        return 0,0
    end
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
    

function census(population)
    dat = broadcast(ind -> [ind.Age ind.G ind.ID], population.individuals)
    return reduce(vcat, dat)
end 
    

function He(population)
    He_acc = 0
    nloci = 0
        
    for ind in population.individuals
            
        for loci in 1:size(ind.genes)[1]
                
            nloci += 1
                
            if all(abs.(ind.genes[loci,:] .- sum(ind.genes[loci,:])/length(ind.genes[loci,:])) .> 10^-2.5 )
                He_acc += 1
            end 
        end 
    end 
        
    return He_acc/nloci
end 

  



########################
## Stochastic slakin ###
########################

mutable struct individualSlakin
    Age::Int64
    G::Float64
    pedigree::Float64
    ID::Float64
    Origin::Int64
    Parents::Tuple{Int64,Int64}
end 

mutable struct stochasticSlakin
    populaitonID::Int64
    individuals::AbstractVector{}
    
    # demographic parameters
    ageStructureParams
    correction 
    p_spawn # probability of spawnign in a given year 
    
    # selection parameters
    s::Float64
    theta::Float64
    gradient::Function
    kernel
    
end 




"""
age_dsn - age distribution
p - frequency at each locus 
"""
function init_individualSlakin(Amax,age_dsn,G_dsn, Origin, pedigree)
    Age = sample(1:Amax,StatsBase.Weights(age_dsn./sum(age_dsn)),1)[1]
 
    G = rand(G_dsn)
    ID = rand()
    parents = (Origin, Origin)
    
    return individualSlakin(Age, G, pedigree, ID, Origin,parents)
end 


function init_stochasticSlakin(ageStructureParams,s,theta, p_spawn, Origin)

    ## compute correction used in deterministic model 
    gradient = x -> exp(-s/2*(x-theta)^2)
    correct = correction(1, s, theta, gradient)

    G_dsn = Distributions.Normal(theta, 1)#V_star(1.0, s))
    
    # get age structure 
    age_dsn = AgeStructuredModels.stable_age_structure(ageStructureParams)    
    
    Neq = floor(Int, sum(age_dsn))  
    individuals = []
    for i in 1:Neq
        
        push!(individuals, init_individualSlakin(ageStructureParams.Amax,age_dsn,G_dsn, Origin, 1.0*Origin))
            
    end 
    kernel = Distributions.Normal(0, sqrt(0.5))
  
    return stochasticSlakin(Origin, individuals, ageStructureParams, correct, p_spawn, s, theta, gradient, kernel)
end



"""
spawning_popualtion - subset of popualtion spawning ins a given year
F - weights given to each indivual inspawning stock 
kernl - distribution genetic effects around midparental value
Origin - Hatcheyr 1 or natrual origin 0
"""
function init_juvinile(spawning_popualtion, F, kernel, Origin)
    # sample parents
    parents = sample(spawning_popualtion, Weights(F./sum(F)), 2, replace=true)

    # sample Haplotype 1   
    G_mean = 0.5*parents[1].G + 0.5*parents[2].G


    # sample trait
    G = G_mean + rand(kernel)
    pedigree = 0.5*parents[1].pedigree + 0.5*parents[2].pedigree
    parents = (parents[1].Origin, parents[2].Origin)
    
    return individualSlakin(1,G,pedigree,rand(),Origin,parents)
end 


"""
Simulates reproduction assuming the denstiy dependent bottle neck occurs befor selection
and immigation of hatchery individuals. 

population - stochasticSlakin
E - environemtnal effect on recruitment 
"""
function reproduction_slakin(population, E)
    x = SSB(population)
    R = population.ageStructureParams.SRCurve(x) * exp(E)
    R *= population.correction
    N = rand(Distributions.Poisson(R))
    
    # get list of mature individuals 
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    mature = F .> 0.0
    spawning_stock = population.individuals[mature]
    F = F[mature]
    # sample a subset that will sapwn this year
    spawning = rand(length(spawning_stock)) .< population.p_spawn
    spawning_stock = spawning_stock[spawning]
    F = F[spawning]
    # init juviniles 
    juviniles = []
    
    if length(spawning_stock) != 0
        for i in 1:N
        
        push!(juviniles, init_juvinile(spawning_stock, F, population.kernel, 0))
            
        end 
        
    end
    
    return juviniles
end 


"""
Simulates reproduction assuming the denstiy dependent bottle neck occurs befor selection
and immigation of hatchery individuals. 

population - stochasticSlakin
E - environemtnal effect on recruitment 
"""
function reproduction_hatchery_slakin(population, Njuviniles, Nspawning)
    
    @assert Nspawning < length(population.individuals)
    
    # hatchery production 
    # get list of mature individuals 
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    mature = F .> 0.0
    spawning_stock = population.individuals[mature]
    F = F[mature]
    # sample a subset that will sapwn this year
    brood_stock_id = sample(1:length(spawning_stock), Nspawning, replace=false)
    brood_stock = spawning_stock[brood_stock_id]
    F =  F[brood_stock_id]
    
    juvinilesH = []
    if length(spawning_stock) != 0
        for i in 1:Njuviniles
        
            push!(juvinilesH, init_juvinile(brood_stock, F, population.kernel, 0))
            
        end 
        
    end
    
    # Natrual production
    x = SSB(population)
    
    R = population.ageStructureParams.SRCurve(x) * exp(E)
    R *= population.correction
    N = rand(Distributions.Poisson(R))
    
    # get list of mature individuals 
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    mature = F .> 0.0
    spawning_stock = population.individuals[mature]
    F = F[mature]
    # sample a subset that will sapwn this year
    spawning = rand(length(spawning_stock)) .< population.p_spawn
    spawning_stock = spawning_stock[spawning]
    F = F[spawning]
    # init juviniles 
    juviniles = []
    
    if length(spawning_stock) != 0
        for i in 1:N
        
        push!(juviniles, init_juvinile(spawning_stock, F, population.kernel, 0))
            
        end 
        
    end
    
    return juvinilesH, juvinlesN
end 



function selection_slakin(juviniles, gradient)
    N = length(juviniles)
    p = broadcast(ind -> gradient(ind.G), juviniles)
    survival = rand(N) .< p
    return juviniles[survival]
end 



function pedigree_dsn(population)
    n = length(population.individuals)
    dsn = zeros(n)
        
    for i in 1:n
        dsn[i] = population.individuals[i].pedigree
    end 
    
    return dsn
end 


function trait_dsn(population)
    n = length(population.individuals)
    dsn = zeros(n)
        
    for i in 1:n
        dsn[i] = population.individuals[i].G
    end 
    
    return dsn
end 

function age_dsn(population)
    n = length(population.individuals)
    dsn = zeros(n)
        
    for i in 1:n
        dsn[i] = population.individuals[i].Age
    end 
    
    return dsn
end 





### single population dynamics
    
function time_step_slakin!(population, E)
    J = reproduction_slakin(population, E)
    J = selection_slakin(J, population.gradient)
    update_population!(population, J)
end 


### single population dynamics
    
function time_step_slakin!(population1, population2, p_sym, p_asym, E1, E2)

        J1 = reproduction_slakin(population1, E1)
        J1 = selection_slakin(J1,population1.gradient)
        
        
        J2 = reproduction_slakin(population2,E2)
        J2 = selection_slakin(J2,population2.gradient)
        
        
        J1, J2 = immigration!(p_sym,J1,J2)
        J1, J2 = directional_immigration!(p_asym,J1, J2)
        
        update_population!(population1, J1)
        update_population!(population2, J2)
end 





end # module 