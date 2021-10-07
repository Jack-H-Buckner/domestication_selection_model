module Normal_population



mutable struct population
    N::Float64 # abundance of popualtion
    r::Float64  
    K::Float64  
    s::Float64
    Vle::Float64 # convolution kernel for Vle 
    order::Int64
    moments::AbstractVector{Float64}
end 




function reproduction!(population, population2)
    Vle = population.Vle
    mu = [population.moments[1],population2.moments[1]]
    sigma = [population.moments[1],population2.moments[1]]
    N = [population.N,population2.N]
    
    mu_prime = sum(mu.*N)/sum(N)
    G = sum(N .* (sigma.^2 .+ mu.^2))/sum(N) - mu_prime^2
    sigma_prime = sqrt(G/2 + Vle/2)
    if population.order == 2
        population.moments = [mu_prime,sigma_prime]
    else
        population.moments  = [mu_prime,sqrt(Vle)]
    end
    population.N = sum(N)
    return population
end 



function selection!(population)
    mu_s = 0
    mu = population.moments[1]
    sigma = population.moments[2]
    N = population.N
    sigma_s = 1/population.s
    
    g = sigma^2*sigma_s^2/(sigma^2 + sigma_s^2)

    m = g*(sigma^2*mu_s + sigma_s^2*mu)/(sigma^2*sigma_s^2)

    N = N*sqrt(g)/(sigma)*exp(-mu^2/(2*sigma^2) - mu_s^2/(2*sigma_s^2) + m^2/(2*g)) #*
    
    population.moments = [m, sqrt(g)]
    population.N = N
    
    return population
end

function density_dependence!(population) 
    r = population.r
    K = population.K
    #N = population.N
    population.N = population.N*exp(r*(1-population.N/K))
    return population
end 



end # module