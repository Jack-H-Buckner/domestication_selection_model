




"""
    AgentBasedSimsUtils

some functions that are useful for defining the write fisher
and stochastic slakin model parameters. I am lumping them off here 
to make the code a bit more readable - this is probabily and excessive 
level of organization.
"""
module AgentBasedSims

# load files and packages 
include("AgeStructuredModels.jl")
include("StockRecruitCurves.jl")
using Distributions
using Roots
using StatsBase


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
    individualWF

State of an individual in the write fisher model and useful statistics.
The state variables are sufficent to define the dynamics of the population 
while the statistics provide useful informaiton for evaluating the outcomes of 
the popualtion dynamics and making statistical inferences. 

State variables: 
ID - unique identifier 
Age - number of time steps since birth 
genes - genetic state of the individual
Origin - hatcheyr or natrual origin 

Statistics: 
G - breading value of individual (summary statistic of genetic state)
pedigree - 
Parents - origin of parents 
"""
mutable struct individualWF # T is an even integer represetning the ploidy of the species of interest
    ID::Float64
    Age::Int64 
    genes::AbstractMatrix{Float64}
    Origin::Int64
    G::Float64
    pedigree::Float64
    Parents::Tuple{Int64,Int64}
end 



"""
    populationWF

State and parameters of the popualton and parameters that govern its dynamics
the state of the population is defined by a list of individuals. In addition
to this list paramters describing the demographic rates of individuals, 
the strength of natrual selection, and the genetic archetecture of individuals 
"""
mutable struct populationWF
    # popualtion stat
    populaitonID::Int64
    individuals::AbstractVector{}
    
    # demographic parameters
    ageStructureParams
    correction 
    p_spawn
    
    # selection parameters
    s::Float64
    theta::Float64
    gradient::Function
    
    # genetic archetecture
    effect::Float64
    Ploidy::Int64 # number of chromospmes
    Nloci::Int64 # number of QTL
end 



"""
    individualSlakin

State of an individual in the Slakin model and useful statistics.
The state variables are sufficent to define the dynamics of the population 
while the statistics provide useful informaiton for evaluating the outcomes of 
the popualtion dynamics and making statistical inferences. 

State variables: 
ID - unique identifier 
Age - number of time steps since birth 
G - breading value of individual
Origin - hatcheyr or natrual origin 

Statistics: 
pedigree - proportion of ancestry from hatcheyr population 
Parents - origin of parents 
"""
mutable struct individualSlakin
    Age::Int64
    G::Float64
    pedigree::Float64
    ID::Float64
    Origin::Int64
    Parents::Tuple{Int64,Int64}
end 


"""
    populationSlakin

State and parameters of the popualton and parameters that govern its dynamics
the state of the population is defined by a list of individuals. In addition
to this list paramters describing the demographic rates of individuals, 
the strength of natrual selection, and the inheritance kernel 
"""
mutable struct populationSlakin
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
    
    # genetic archetecture 
    kernel
    
end 


"""
     init_individual(General params,model specific params)

Initalizes an individual for the agent based model. Uses multiple dispatch 
to initialize individuals for either the Slakin or WF model.

General params:
Amax - Maximum age 
age_dsn - Age distribution
Origin - populaiton tag (0 or 1)

WF params:
p - allel frequency at each locus 
E - Effect size of each allel
Ploidy - number of chromosomes
Nloci - number of QTLs

Slakin params:
G_dsn - distribution of breading values 

"""
function init_individual(Amax,age_dsn,Origin, p,E,Ploidy,Nloci)
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
  
    return individualWF(rand(),Age,genes,Origin,G,Origin,(Origin,Origin))
        
end 
function init_individual(Amax,age_dsn,Origin, G_dsn)
    Age = sample(1:Amax,StatsBase.Weights(age_dsn./sum(age_dsn)),1)[1]
 
    G = rand(G_dsn)
    ID = rand()
    parents = (Origin, Origin)
    
    return individualSlakin(Age, G, 1.0*Origin, ID, Origin,parents)
end 


"""
    init_population(general params, model specific params)

Initializes a populaiton object and used multiple dispatch to 
choose between Slakin and Write Fisher model. 

General parameters:
ageStructureParams - parameters that define demographic rates
p_spawn - probabilty an individuals spawns in a given year
Origin - population of origin (0 or 1)
s - selection strength
theta - optimal trait value
 

Write Fisher parameters:
Ploidy - number of chromosomes
Nloci - number of QTLs

Slakin parameters:
no additonal params inheritance kernels set to a default of 
Normal(0,0.5)
"""
function init_population(ageStructureParams,p_spawn,Origin,s,theta,Ploidy,Nloci)
    
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
    
        push!(individuals, init_individual(ageStructureParams.Amax,age_dsn, Origin, p,E,Ploidy,Nloci))
            
    end 
        
    return populationWF(Origin,individuals, ageStructureParams, correct, p_spawn, s,theta,gradient,E,Ploidy,Nloci)
    
end    
function init_population(ageStructureParams,p_spawn,Origin,s,theta)

    ## compute correction used in deterministic model 
    gradient = x -> exp(-s/2*(x-theta)^2)
    correct = correction(1, s, theta, gradient)

    G_dsn = Distributions.Normal(theta, 1)#V_star(1.0, s))
    
    # get age structure 
    age_dsn = AgeStructuredModels.stable_age_structure(ageStructureParams)    
    
    Neq = floor(Int, sum(age_dsn))  
    individuals = []
    for i in 1:Neq
        
        push!(individuals, init_individual(ageStructureParams.Amax,age_dsn,Origin,G_dsn))
            
    end 
    kernel = Distributions.Normal(0, sqrt(0.5))
  
    return populationSlakin(Origin, individuals, ageStructureParams, correct, p_spawn, s, theta, gradient, kernel)
    
end





"""
    init_juvinile(spawning_stock,F, Origin)

Samples parents at random from the spawning stock and initializes a
age one individual to add to the population. Multiple dispatch is used
to determin if Slakin of WF model is applied. 

the inheritance kernel is supplied as an addtional prameter for the Slakin model 
"""
function init_juvinile(spawning_stock,F, Origin,Ploidy::Int64)
    # sample parents
    parents = sample(spawning_stock, Weights(F./sum(F)), 2, replace=true)
    
    # genetics 
    H1 = mapslices(x -> rand(x,floor(Int,Ploidy/2)), parents[1].genes, dims = 2)
    H2 = mapslices(x -> rand(x,floor(Int,Ploidy/2)), parents[2].genes, dims = 2)
    genes = hcat(H1,H2)
    G = sum(genes)
    
    # pedigree
    pedigree = 0.5*parents[1].pedigree + 0.5*parents[2].pedigree
    parents = (parents[1].Origin, parents[2].Origin)
    
    return individualWF(rand(),1, genes, Origin, G, pedigree,  parents)
end 



function init_juvinile(spawning_stock,F, Origin, kernel)
    # sample parents
    parents = sample(spawning_stock, Weights(F./sum(F)), 2, replace=true)
    
    # genetics 
    G = 0.5*parents[1].G + 0.5*parents[2].G + rand(kernel)
    
    # pedigree
    pedigree = 0.5*parents[1].pedigree + 0.5*parents[2].pedigree
    parents = (parents[1].Origin, parents[2].Origin)
    
    return individualSlakin(1,  G, pedigree, rand(1)[1], Origin, parents)
end 

"""
    init_juvinile_hatchery(N, G_dsn)
"""
function init_juvinile_hatchery(N, G_dsn)
    
    juviniles = []
    
    for i in 1:N
        J = individualSlakin(1,rand(G_dsn), 1.0, rand(),1, (1,1))
        push!(juviniles, J)
    end 
    
    return juviniles
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


"""
    hatchery_removals(population, E)

removes hathcery fish from popualtion with probability exp(-E)
"""
function hatchery_removals!(population, E)
    
    N = length(population.individuals)
    survival = ones(N)
        
    for i in 1:N
        if population.individuals[i].Origin == 1
            if exp(-E) < rand()
                survival[i] = 0
            end 
        end
    end 
    population.individuals = population.individuals[floor.(Int,survival)]
end 




"""
Simulates reproduction assuming the denstiy dependent bottle neck occurs befor selection
and immigation of hatchery individuals. 

population - stochasticSlakin
E - environemtnal effect on recruitment 
"""
function reproduction(population::populationSlakin, E)
    x = SSB(population)
    R = population.ageStructureParams.SRCurve(x) * exp(E)
    R *= population.correction
    N = rand(Distributions.Poisson(R))
    
    # get list of mature individuals 
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    mature = F .> 0.0
    spawning_stock =population.individuals[mature]
    F = F[mature]
    
    # sample a subset that will sapwn this year
    spawning = rand(length(spawning_stock)) .< population.p_spawn
    spawning_stock = spawning_stock[spawning]
    F = F[spawning]
    
    # init juviniles 
    juviniles = []
    
    if length(spawning_stock) != 0
        for i in 1:N
        
            push!(juviniles, init_juvinile(spawning_stock,F, population.populaitonID, population.kernel))
            
        end 
        
    end
    
    return juviniles
end 
function reproduction(population::populationWF, E)
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
        
            push!(juviniles, init_juvinile(spawning_stock,F, population.populaitonID, population.Ploidy))
            
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
function reproduction_hatchery(population::populationSlakin, Njuviniles, Nspawning, E)
    
    if 0.5*Nspawning > length(population.individuals)
        
        Nspawning = floor(Int, Nspawning/2)
        
    end
    
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
        
            push!(juvinilesH, init_juvinile(brood_stock,F, population.populaitonID, population.kernel))
                
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
    juvinilesN = []
    
    if length(spawning_stock) != 0
        for i in 1:N
        
            push!(juvinilesN, init_juvinile(spawning_stock,F, 1, population.kernel))
            
        end 
        
    end
    
    return juvinilesH, juvinilesN
end 
function reproduction_hatchery(population::populationWF, Njuviniles, Nspawning, E)
    
    if 0.5*Nspawning > length(population.individuals)
        
        Nspawning = floor(Int, Nspawning/2)
        
    end
    
    # hatchery production 
    # get list of mature individuals 
    F = broadcast(ind -> population.ageStructureParams.Fecundity[ind.Age],  population.individuals)
    mature = F .> 0.0
    spawning_stock = population.individuals[mature]
    F = F[mature]
    
    # sample a subset that will spawn this year
    brood_stock_id = sample(1:length(spawning_stock), Nspawning, replace=false)
    brood_stock = spawning_stock[brood_stock_id]
    F =  F[brood_stock_id]
    
    juvinilesH = []
    if length(spawning_stock) != 0
        for i in 1:Njuviniles
            
            push!(juvinilesH, init_juvinile(brood_stock,F, population.populaitonID, population.Ploidy))
            
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
    juvinilesN = []
    
    if length(spawning_stock) != 0
        for i in 1:N
        
            push!(juvinilesN, init_juvinile(spawning_stock, F, 0, population.Ploidy))
            
        end 
        
    end
    
    return juvinilesH, juvinilesN
end 


function selection(gradient, juviniles)
    survival = zeros(length(juviniles))
    i=0
    for ind in juviniles
        i+=1
        if gradient(ind.G) > rand(1)[1]
            survival[i] = 1
        else
            survival[i] = 0
        end
    end 
    return juviniles[survival.==1]
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

end 
    
    
### single population dynamics

"""
    time_step!(population, optional params)

Updates the population  
"""
function time_step!(population, E)

    J = reproduction(population, E)
    J = selection(population.gradient, J)
    update_population!(population, J)

end 

# with immigration  
function time_step!(population, Nim, Gdsn, E, before)
    if before
        J = reproduction(population, E )
        J_im = init_juvinile_hatchery(Nim, Gdsn)
        J = selection(population.gradient, vcat(J,J_im))
        update_population!(population, J)
    else
        J = reproduction(population, E )
        J = selection(population.gradient, J)
        J_im = init_juvinile_hatchery(Nim, Gdsn)
        update_population!(population, vcat(J,J_im))  
    end
end 


# with hatchery reproduction   
function time_step!(population, Hgrad, Njuv, Nspawn, E)
    JH,JN = reproduction_hatchery(population, Njuviniles, Nspawning, E)
    JN = selection(population.gradient, JN)
    JH = selection(Hgrad, JN)
    update_population!(population, vcat(JN,JH))    
end 



# two population dynamics 
function time_step!(population1, population2, before, Njuv, Nspawn, p_sym, p_asym, E1, E2)
    
    if before
        J1 = reproduction(population1, E1)
        
        J2H,J2N = reproduction_hatchery(population, Njuviniles, Nspawning, E)

        J1, J2N = immigration!(p_sym,J1,J2N)
        J1, J2H = directional_immigration!(p_asym,J1, J2H)

        J1 = selection(population1.gradient,J1)
        J2H = selection(population2.gradient,J2N)
        J2N = selection(population2.gradient,J2H)
        
        
        update_population!(population1, J1)
        update_population!(population2, vcat(J2N,J2H))
    else
        J1 = reproduction(population1, E1)
        J1 = selection(population1.gradient,J1)


        J2H,J2N = reproduction_hatchery(population, Njuviniles, Nspawning, E)
        J2H = selection(population2.gradient,J2N)
        J2N = selection(population2.gradient,J2H)


        J1, J2N = immigration!(p_sym,J1,J2N)
        J1, J2H = directional_immigration!(p_asym,J1, J2H)

        update_population!(population1, J1)
        update_population!(population2, vcat(J2N,J2H))
            
    end
end 





## population statistics

"""
    trait_dsn(population)

distribution of breeding values in the population 
"""
function trait_dsn(population)
    Gls = []
    for ind in population.individuals
        push!(Gls, ind.G)
    end 
    return Gls 
end 

"""
    trait_dsn_moments(population)

mean and variance of breeding value distribution 
"""
function trait_dsn_moments(population)
    
    Gls = trait_distribution(population)
    if length(Gls) > 0 
        mu = sum(Gls)/length(Gls)
        sigma = sqrt(sum((Gls .- mu).^2 ./(length(Gls) - 1)))
        return mu, sigma
    else
        return 0,0
    end
end 


"""
    trait_dsn_age1(population)

distribution of breeding values among recruits
"""
function trait_dsn_age1(population)
    Gls = []
    for ind in population.individuals
        if ind.Age == 1
            push!(Gls, ind.G)
        end
    end 
    return Gls 
end 


"""
    trait_dsn_age1(population)

Mean and variance of breeding values among recruits
"""
function trait_dsn_moments_age1(population)
    Gls = trait_distribution_age1(population)
    mu = sum(Gls)/length(Gls)
    sigma = sqrt(sum((Gls .- mu).^2 ./(length(Gls) - 1)))
    return mu, sigma
end 



"""
    pedigree_dsn(population)

list of state varibles for each individual in the population 
"""
function census(population)
    dat = broadcast(ind -> [ind.Age ind.G ind.ID ind.pedigree ind.Origin], population.individuals)
    return reduce(vcat, dat)
end 
    


"""
    pedigree_dsn(population)

Distribution of pedigree values in the population 
"""
function pedigree_dsn(population)
    n = length(population.individuals)
    dsn = zeros(n)
        
    for i in 1:n
        dsn[i] = population.individuals[i].pedigree
    end 
    
    return dsn
end 



"""
    age_dsn(population)

Returns the age distribution of the population 
"""
function age_dsn(population)
    n = length(population.individuals)
    dsn = zeros(n)
        
    for i in 1:n
        dsn[i] = population.individuals[i].Age
    end 
    
    return dsn
end 



"""
    He(population)

average heterozygosity of a locus in the population 
"""
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


end # module 


