module i_state_dsn_pedigree_trait

include("utils.jl")
include("../age_structure/age_structure_model.jl")
using Roots
using Distributions
using FFTW
using DSP
using Plots

mutable struct population
    abundance::AbstractVector{Float64}
    trait::AbstractVector{AbstractMatrix{Float64}} # columns are ages 
    grid_z::AbstractVector{Float64} # nodes where trait is defined
    grid_p::AbstractVector{Float64} # nodes where trait is defined
    fecundity::AbstractVector{Float64}
    survival::AbstractVector{Float64}
    A_max::Int64 # max age
    a::Float64 # stock recruit curve - based on total fecundity
    b::Float64
    gradient::AbstractVector{Float64} # age dependent seleciton gradient columsn are ages
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel 
    n::Int64
end 

function init_population(A_max, survival, fecundity, r , K, theta, s, min, max, dx, n, Vle)
    
    # set age distribution - equilibrium 
    LEP = age_structure_model.LEP(1.0, survival, fecundity, A_max)
    
    a = r/LEP
    b = a/K
    f(x) = a*x/(1+b*x)
    g(x )= f(x) - x/LEP

    x_eq = Roots.find_zero(g, (10^-6, 10^6*K*LEP))
    y_eq = f(x_eq)
    
    abundance = zeros(A_max)
    N = y_eq
    for i in 1:A_max
        abundance[i] = N
        N *= survival[i]
    end 
  
    # set trait at optimum value and variance = Vle
    grid_z = collect(min:dx:max)
    grid_p = collect(1:(2^n+1))
  
    d = Distributions.Normal(theta, sqrt(Vle))
    trait = repeat([zeros(length(grid_p),length(grid_z))],A_max)
    
    for i in 1:A_max
        trait[i][1,:] = pdf.(d, grid_z) ./ sum(pdf.(d, grid_z))
    end 

    
    
    # Vle
    de = Distributions.Normal(theta, sqrt(Vle/2))
    grid_Vle = collect((theta-3*Vle):dx:(theta+3*Vle))
    Vle = pdf.(de,grid_Vle)
    m = length(Vle)
    
    # gradient
    gradient = exp.(-s/2*(grid_z .- theta).^2)
    
    pop = population(abundance, trait, grid_z, grid_p, fecundity, survival,A_max, a, b, gradient,m,Vle, 2^n+1)
    
    return pop
end 



"""
returns trait distribution for new age class and 
the total spawning stock fecundity.
"""
function reproduction(population)

    f_total = sum(population.abundance .* population.fecundity)
    dsn = population.trait .* (population.abundance .* population.fecundity) ./ f_total
    dsn = sum(population.trait)
    dsn = dsn ./ sum(dsn)

    s = size(dsn)
    dsn_zs = zeros(2*s[1], 2*s[2])
    dsn_zs[1:s[1],1:s[2]] = dsn

    dsn_zs = real(ifft(fft(dsn_zs).*fft(dsn_zs)))
    dsn = utils.colapse_M(dsn_zs)
    dsn = dsn./sum(dsn)
    


    N = size(dsn)[2]
    m = convert(Int64,floor(population.m/2+1))
    dsn = transpose(DSP.conv(transpose(dsn), population.Vle))
     

    dsn= dsn[:,m:(N+m-1)]
    dsn = dsn ./ sum(dsn)

 
    return dsn, f_total
end 




"""
selection on juviniles 
"""
function selection(dsn, N, population)
    dsn = transpose(transpose(dsn) .* population.gradient)
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
    return R
end 



"""
Immigration of juviniles
"""
function immigration(dsn, N, dsn_im, N_im, populations)
    N_total = N + N_im
    p = N/N_total
        
    dsn = p.*dsn 
    dsn[populations.n,:] = (1-p).* dsn_im
    return dsn, N_total
end 
    
"""
Updates the age structure of the popuatlion and adds recruits
"""
function ageing!(population, R, dsn_R)
    A_max = population.A_max
    population.abundance = population.abundance .* population.survival
    plus_group = sum(population.abundance[end-1:end])
    new_A = zeros(A_max)
    new_A[1] = R 
    new_A[end] = plus_group
    

    
    new_A[2:end-1] = population.abundance[1:end-2]
    
    s = size(population.trait[A_max])
    new_dsn = repeat([zeros(s[1],s[2])],A_max)
    tplus = copy(population.trait[A_max])
    tplus_m1 = copy(population.trait[A_max-1])
    
    p = population.abundance[end]./plus_group
    tplus = p*tplus .+ (1-p)*tplus_m1
    new_dsn[A_max] = tplus
        
    new_dsn[2:end-1] = population.trait[1:end-2]
        
    new_dsn[1] = dsn_R
    population.abundance = new_A
    population.trait = new_dsn
    
    return population
end 


    
end # module 