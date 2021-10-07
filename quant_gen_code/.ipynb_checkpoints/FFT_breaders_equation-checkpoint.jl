module FFT_breaders_equation

using FFTW

mutable struct population
    N::Int64 # abundance of popualtion 
    m::Int64 # length of convolution kernel 
    Vle::AbstractVector{Float64} # convolution kernel for Vle 
    values::AbstractVector{Float64} # trait distribution 
    grid::AbstractVector{Float64} # grid values for trait distribution
end 

function reproduction!(population)
    dsn = population.values
    v1 = DSP.conv(dsn, dsn)[1:2:(2*N+1)]
    v1 = v1./sum(v1)
    dsn_e = population.values[Vle]
    v1 = DSP.conv(v1, dsn_e)[(m+1):(N+m+1)]
    population.values = v1
    return population
end 

end # module