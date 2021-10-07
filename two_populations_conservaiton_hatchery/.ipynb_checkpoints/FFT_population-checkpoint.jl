module FFT_population

using FFTW
using DSP
using Distributions
mutable struct population
    N::Float64 # abundance of popualtion
    r::Float64  
    K::Float64  
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel for Vle 
    values::AbstractVector{Float64} # trait distribution 
    gradient::AbstractVector{Float64} # survival probility under selection
    grid::AbstractVector{Float64} # grid values for trait distribution
end 

function reproduction!(population)
    dsn = population.values
    N = length(population.grid)-1
    v1 = DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    v1 = v1./sum(v1)
    dsn_e = population.Vle
    m = population.m
    v1 = DSP.conv(v1, dsn_e)[(m+1):(N+m+1)]
    population.values = v1 ./ sum(v1)
    
    return population
end 


function flow!(population, dsn, N)
    N_old = population.N
    population.N += N
    population.values = (N_old*population.values .+ N*dsn)./population.N
    return population
end 


function selection!(population)
    population.values .*= population.gradient
    population.N *= sum(population.values)
    population.values = population.values ./ sum(population.values)
    return population
end

function density_dependence!(population) 
    r = population.r
    K = population.K
    #N = population.N
    population.N = population.N*exp(r*(1-population.N/K))
    return population
end 


function init_pop(trait_diff,r,K,s)

    N = 500
    N0 = 1.0
    x_e = -1.0:(20/N):1.0
    Vle_val = 1
    Vle = pdf.(Distributions.Normal(0,Vle_val), x_e)
    Vle = Vle./sum(Vle)
    N_e = length(x_e)
    m = convert(Int64,floor(N_e/2))
    grid = -10.0:(30/N):20.0
    values = pdf.(Distributions.Normal(0,Vle_val),grid)
    values = values./sum(values)
    gradient = exp.(- s/2*grid.^2)
    dsn_flow =  pdf.(Distributions.Normal(trait_diff ,Vle_val),grid)
    dsn_flow = dsn_flow  ./ sum(dsn_flow )
    return population(N0,r,K,m,Vle,values,gradient,grid), dsn_flow
end 





end # module