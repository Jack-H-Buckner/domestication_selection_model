module AgeTraitIntraCohortComp

include("AgeStructuredModels.jl")
include("AgeTraitStructuredModels.jl")
using Distributions
using DSP

mutable struct params
    k::Float64
    Rstar::Float64
    survival::AbstractVector{Float64}
    fecundity::AbstractVector{Float64}
    r::Float64
    Wa::AbstractVector{Float64}
    functional_form::String
    Vle::Float64
    theta::Float64
    s::Float64
    min::Float64
    max::Float64
    dx::Float64
end 


mutable struct population
    # states
    abundanceN::AbstractVector{Float64}
    abundanceH::AbstractVector{Float64}
    traitN::AbstractMatrix{Float64} # columns are ages 
    traitH::AbstractMatrix{Float64}
    grid::AbstractVector{Float64} # nodes where trait is defined
    
    # dempgraphic paramters
    survival::AbstractVector{Float64}
    fecundity::AbstractVector{Float64}
    functional_form::String
    b0::Float64
    ba::AbstractVector{Float64}
    k
    Rstar
    LEP
    # genetic parameters
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
    correction::Float64
    
    # if selection is before desntiy dependnece then use corection to get to k*Rstar at equilibrium
    # if seleciton is after denstiy dependnece the use correction to get to Rstar at equilibirum
    
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



function init(params)
    k = params.k; Rstar = params.Rstar; survival = params.survival; fecundity = params.fecundity
    r = params.r;  Wa = params.Wa; functional_form = params.functional_form; Vle = params.Vle
    theta = params.theta; s = params.s; min = params.min; max = params.max; dx = params.dx
    Amax = length(survival)
    #### set up genetic paramters and initialize genetic states ####
    
    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(AgeTraitStructuredModels.V_star(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    traitN = transpose(repeat(transpose(trait),Amax ))
    
    traitH = pdf.(d, grid)
    traitH = traitH ./sum(traitH)
    traitH = transpose(repeat(transpose(traitH),Amax ))
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    Vle_ = Vle
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    d = Distributions.Normal(theta, sqrt(AgeTraitStructuredModels.V_star_prime(Vle_, s)))
    trait = pdf.(d, grid)
    trait = trait ./ sum(trait)
    gradient = exp.(-s/2*(grid .- theta).^2)
    correction = 1/sum(trait .* gradient)
    

    #### initialize demographic paramters ####
    
    # base paramters
    Sa = vcat([1], broadcast(i -> prod(survival[1:i]), 1:(Amax-1)))
    b0 = Wa[1]
    ba = vcat(Wa[2:end], [0]) .* (1-r).^vcat(Wa[2:end], [0])
    LEP = sum(Sa.*fecundity)
    # calcualte scale
    if functional_form =="Ricker"
        b = log(k)/(b0*k*Rstar+sum(Rstar*Sa.*ba))
    elseif functional_form == "BevertonHolt"
        b = (k-1)/(b0*k*Rstar+sum(Rstar*Sa.*ba))
    end 
    b0 *= b
    ba .*= b
    #### initialize demographic state ####

    abundanceN = Rstar*Sa
    abundanceH = zeros(Amax)
    
    return population(  abundanceN,abundanceH,
                        traitN,traitH, grid,
                        survival, fecundity,
                        functional_form, b0, ba,k,Rstar,LEP,
                        gradient, m, Vle,
                        correction,
                        Vle_, s, theta)
    
end 



function init_imigrants(population, N, mean)
    grid = population.grid
    d = Distributions.Normal(mean, sqrt(population.Vle_))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    return immigrants(N,trait,grid)
end 




function reproduction(population)
    # useful paramters
    k = population.k
    Rstar = population.Rstar
  
    
    ### egg production ###
    f_totalN = sum(population.abundanceN .* population.fecundity)
    f_totalH = sum(population.abundanceH .* population.fecundity)
    
    
    ### new distribution ###
    
    # contributions
    dsnN = population.traitN * (population.abundanceN .*  population.fecundity)# ./ f_totalN
    dsnH = population.traitH * (population.abundanceH .*  population.fecundity) #./ f_totalH
    dsn = (dsnN + dsnH)./ (f_totalN+f_totalH)
    
    
    # convolution
    N = length(population.grid)-1
    dsn = DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    dsn = dsn./sum(dsn)
    
    # segregation variance
    m = convert(Int64,floor(population.m/2+1))
    dsn = DSP.conv(dsn, population.Vle)
    dsn= dsn[m:(N+m)]
    
    # rescale
    dsn = dsn ./ sum(dsn)
    
    
    return dsn, k*(f_totalN + f_totalH)/population.LEP
end 



function selection(dsn, N, population)
    dsn = dsn .* population.gradient
    survival = sum(dsn)
    N = N * survival
    dsn = dsn ./ survival
    return dsn, N
end 



function density_dependence(N, H, population)
    b0 = population.b0
    ba = population.ba
    
    J = N + H
    Na = population.abundanceN .+ population.abundanceH
    if population.functional_form == "Ricker"
        p_survival = exp(-b0*J-sum(Na.*ba))
    else
        p_survival = 1/(1 + b0*J + sum(Na.*ba))
    end 
    return p_survival
end 


function ageing!(population, R, dsn_R, NH, dsnH)
    
    population.abundanceN = population.abundanceN .* population.survival
    population.abundanceH = population.abundanceH .* population.survival
    
    Amax = length(population.survival) 
    
    new_N = zeros(Amax )
    new_N[1] = R 
    new_N[2:end] = population.abundanceN[1:end-1]
    
    new_H = zeros(Amax )
    new_H[1] = NH
    new_H[2:end] = population.abundanceH[1:end-1]
    
    
    N = length(dsn_R)
    new_dsnN = zeros(N,Amax )
    new_dsnN[:,2:end] = population.traitN[:,1:end-1]
    new_dsnN[:,1] = dsn_R
    
    new_dsnH = zeros(N,Amax )
    new_dsnH[:,2:end] = population.traitH[:,1:end-1]
    new_dsnH[:,1] = dsnH

    
    population.abundanceN = new_N
    population.traitN = new_dsnN
    
    population.abundanceH = new_H
    population.traitH = new_dsnH
end 





function time_step_DSI!(population, immigrants)
    # reproduction 
    dsnN, N = reproduction(population)
    
    # densit y dependence 
    N *= density_dependence(N, 0, population)
    
    # selection
    dsnN, N = selection(dsnN, N, population)
    N *= population.correction
    
    # immigration 
    H = immigrants.N
    dsnH = immigrants.trait
    
    # aging
    ageing!(population, N, dsnN, H, dsnH)
    
end 




function time_step_IDS!(population, immigrants)
    # reproduction 
    dsnN, N = reproduction(population)
    
    # immigration 
    H = immigrants.N
    dsnH = immigrants.dsn
    
    # density dependence 
    N *= density_dependence(N, H, population)
    
    # selection
    dsnN, N = selection(dsnN, N, population)
    N *= population.correction
    
    dsnH, H = selection(dsnH, H, population)
    H *= population.correction
    
    
    # aging
    ageing!(population, N, dsnN, H, dsnH)
    
end 


# statistics
function fittness(population)
    dsn, f = reproduction(population)
    W = sum(dsn .* population.gradient)
    return W
end 



end 