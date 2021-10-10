"""
ths file provides funciton that define an age structured popualtion model 
that tracks he full age dependent distribution of a trait in the populaiton. 

The code is built around a centeral mutable struct, 'population'. Popualtion stores 
a vector with abundcaes of each age class, and an array that stores the density of the 
trait distribution for each cohort at a set of nodes. In additon, to these states the object
stores paramters that describe the fecudndity and survival of each age class, a stock - recruit
relationship, selection gradient and trait variance. 

Around this basic structure, several methods are defined to update the state of the popuatlion
as as differnt processes occur, such as immigraiton, reproduction, aging and recruitment. 
"""
module full_trait_distribution

using Distributions
using DSP
using Plots


mutable struct population
    abundance::AbstractVector{Float64}
    trait::AbstractMatrix{Float64} # columns are ages 
    grid::AbstractVector{Float64} # nodes where trait is defined
    fecundity::AbstractVector{Float64}
    survival::AbstractVector{Float64}
    A_max::Int64 # max age
    r::Float64 # stock recruit curve - based on total fecundity
    K::Float64
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
end 

"""
returns trait distribution for new age class and 
the total spawning stock fecundity.
"""
function reproduction(population)

    f_total = sum(population.abundance .* population.fecundity)
    dsn = population.trait * (population.abundance .* population.fecundity) ./ f_total
    

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
"""
function recruitment(f_total, population )
    R = population.r*f_total/(1+population.r*f_total/population.K )
    return R
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
function ageing!(populaiton, R, dsn_R)
    population.abundance = population.abundance .* population.survival
    plus_group = sum(population.abundance[end-1:end])
    new_A = zeros(population.A_max + 1)
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


function init_population(A_max, survival, fecundity, r ,K, theta, s, min, max, dx, Vle)
    
    # set age distribution - ad hoc method
    abundance = zeros(A_max)
    N = 1
    for i in 1:A_max
        abundance[i] = N
        N *= survival[i]
    end 
    F = fecundity.*abundance
    N0 = F*K*(r+1.0)./r
    abundance .*= N0
    
    # set trait at optimum value and variance = Vle
    grid = collect(min:dx:max)
    d = Distributions.Normal(theta, sqrt(Vle))
    trait = pdf.(d, grid)
    trait = trait ./sum(trait)
    trait = transpose(repeat(transpose(trait),A_max))
    
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    gradient = exp.(-s/2*(grid .- theta).^2)
    
    pop = population(abundance, trait, grid, fecundity, survival,A_max, r, K, gradient,m,Vle)
    
    return pop
end 

end # module 