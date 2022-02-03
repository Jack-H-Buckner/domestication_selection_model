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
module FullTriatDistributon

include("age_structure_model.jl")
using Distributions
using DSP
using Plots
using Roots
using FFTW

mutable struct population
    abundance::AbstractVector{Float64}
    trait::AbstractMatrix{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
    fecundity::AbstractVector{Float64}
    survival::AbstractVector{Float64}
    A_max::Int64 # max age
    a::Float64 # stock recruit curve - based on total fecundity
    b::Float64
    correction::Float64
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
end 

"""
returns trait distribution for new age class and 
the total spawning stock fecundity.
"""
function reproduction_fft(population)

    f_total = sum(population.abundance .* population.fecundity)
    dsn = population.trait * (population.abundance .* population.fecundity) ./ f_total
    

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

    f_total = sum(population.abundance .* population.fecundity)
    dsn = population.trait * (population.abundance .* population.fecundity) ./ f_total
    

    #Plots.plot!(population.grid,dsn )
    # convolution - random mating 
    N = length(population.grid)-1
    dsn_zs = vcat(dsn,zeros(length(dsn)))
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
"""
function recruitment(f_total, population )
    R = population.a*f_total/(1+population.b*f_total )
    return R .* population.correction
end 



"""
Immigration of juviniles
"""
function immigration(dsn, N, dsn_im, N_im)
    N_total = N + N_im
    p = N/N_total
    dsn = p.*dsn .+ (1-p).* dsn_im
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
    population.abundance = population.abundance .* population.survival
    plus_group = sum(population.abundance[end-1:end])
    new_A = zeros(population.A_max)
    new_A[1] = R 
    new_A[end] = plus_group
    

    
    new_A[2:end-1] = population.abundance[1:end-2]
    
    N = length(population.grid)
    new_dsn = zeros(N,population.A_max)
    new_dsn[:,1] = dsn_R
    plus_group_dsn = (population.abundance[end-1]*population.trait[:,end-1] .+ population.abundance[end]*population.trait[:,end]) ./plus_group
    new_dsn[:,end] = plus_group_dsn
    
    new_dsn[:,2:end-1] = population.trait[:,1:end-2]
    
    population.abundance = new_A
    population.trait = new_dsn
    
    return population
end 


function V_star(Vle, s)
    sigma_s = 1/s
    V_prime = V -> (V*sigma_s^2/(V + sigma_s^2))/2 + Vle/2
    V0 = Vle
    V1 = V_prime(Vle)
    while abs(V0-V1) > 10^-6
        V0 = V1
        V1 = V_prime(Vle)
    end 
    return V1
end 

function init_population(A_max, survival, fecundity, r , K, theta, s, min, max, dx, Vle)
    
    # set age distribution - equilibrium 
    LEP = age_structure_model.LEP(1.0, survival, fecundity, A_max)
    
    a = r/LEP
    b = a/K
    f(x) = a*x/(1+b*x)
    g(x )= f(x) - x/LEP
    print(LEP)
    x_eq = Roots.find_zero(g, (10^-6, 10^6*K*LEP))
    y_eq = f(x_eq)
    
    abundance = zeros(A_max)
    N = y_eq
    for i in 1:A_max
        abundance[i] = N
        N *= survival[i]
    end 

    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),A_max))
    
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    gradient = exp.(-s/2*(grid .- theta).^2)
    correction = 1/sum(trait .* gradient)
    
    pop = population(abundance, trait_A, grid, fecundity, survival,A_max, a, b, correction, gradient,m,Vle)
    
    return pop
end 



function init_population_2(A_max, survival, fecundity, R_star, alpha , beta, theta, s, min, max, dx, Vle, correction)
    
    # set age distribution - equilibrium 
    LEP = age_structure_model.LEP(1.0, survival, fecundity, A_max)
    
    
    abundance = zeros(A_max)
    N = R_star
    for i in 1:A_max
        abundance[i] = N 
        N *= survival[i]
    end 
    # set trait at optimum value and variance = Vle 
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(V_star(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    trait_A = transpose(repeat(transpose(trait),A_max))
    
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    if correction
        gradient = exp.(-s/2*(grid .- theta).^2)
        correction = 1/sum(trait .* gradient)
    else
        correction = 1
    end 
    
    
    
    pop = population(abundance, trait_A, grid, fecundity, survival,A_max, alpha, beta, correction, gradient,m,Vle)

    return pop
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

end # module 