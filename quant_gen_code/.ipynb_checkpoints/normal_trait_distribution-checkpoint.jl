module normal_trait_distribution
using Distributions

"""
reproduction(x, mu, sigma, f, V_le)
Models the  mean trait values produced by random mating 
of a group of sub popuatlions.

Arguments:
mu - vector of mean trait values from each contributing populaiton
sigma - vector of trait variances from each contributing populaiton
x - vector of sizes of each contribution population
f - float64 fecundity
V_ls - float64 variance at likage equilibrium

Value:
x_prime - float64 density of new juviniles
mu_prime - float64 - mean trait value 
sigma_prime - float64 - variance in trait values
"""
function reproduction(x, mu, sigma, f, V_le)
    @assert length(mu) == length(sigma)
    @assert length(x) == length(mu)
    mu_prime = sum(mu.*x)/sum(x)
    G = sum(x .* (sigma.^2 .+ mu.^2))/sum(x) - mu_prime^2
    sigma_prime = G/2 + V_le/2
    x_prime = f*sum(x)
    if sigma_prime < 0
        print("issue")
        print("\n")
        sigma_prime = 0.001
    end
    return x_prime, mu_prime, sqrt(sigma_prime) 
end 

"""
selection(x, mu, sigma, mu_s, sigma_s)
Models the effect of selection on the trait distribution of a popualtion

Arguments:
x - float64 populaiton density
mu - float64 mean trait value 
sigma - float64 variance of trait
mu_s - optimal trait value
sigma_s - variance of seleciton 

Value:
x_prime - float64 size of popualtion after selection
mu_prime - float64 mean triat value after selection
sigma_prime - float64 variacne after selection
"""
function selection(x, mu, sigma, mu_s, sigma_s)
    
    g = sigma^2*sigma_s^2/(sigma^2 + sigma_s^2)
    
    m = g*(sigma^2*mu_s + sigma_s^2*mu)/(sigma^2*sigma_s^2)
    
    x = x*sqrt(g)/(sigma)*exp(-mu^2/(2*sigma^2) - mu_s^2/(2*sigma_s^2) + m^2/(2*g)) #*
    return x, m, sqrt(g)
end 


function selection_vector(x, mu, sigma, mu_s, sigma_s)
    g = sigma^2*sigma_s^2/(sigma^2 + sigma_s^2)
    m = g*(sigma^2*mu_s + sigma_s^2*mu)/(sigma^2*sigma_s^2)
    x = x*sqrt(g)/(sigma)*exp(-mu^2/(2*sigma^2) - mu_s^2/(2*sigma_s^2) + m^2/(2*g)) #*
    return [x, m, sqrt(g)]
end 

"""
competition(x,k)
Models density dependent mortality

Arguments:
x - vector popualtion sizes
k - strength of denstiy dependence

Value:
x - vector with popualtion sizes of each populaiton after density dependence
"""
function competition(x,f,k)
    h = sum(x)*(f-1)/(k*f)
    return x./(1+h)
end 


"""
brood_stock(x,max)
calcualtes size of brood stock give popualtion size x

Arguments:
x - float64 popualation size
prop - float64 proportion ofpopualtiont taken ad brood stock
max - maxumum number of individuals taken as brood stock

Value:
x - size of population after brood stock is taken
B - number of individuals taken as brood stock
"""
function brood_stock(x,prop,max)
    if prop*x < max
        return (1-prop)*x, prop*x
    else
        return x - max, max
    end
end 


"""
reproduction(x, mu, sigma, f, V_le)
Models the  mean trait values produced by random mating 
of a group of sub popuatlions. Each one is assumed to be finite. 

The mean and variacne of the trait vlues are updated before reproduction
to reflect stochasticiy in the survival process. 

Arguments:
mu - vector of mean trait values from each contributing populaiton
sigma - vector of trait variances from each contributing populaiton
x - vector of sizes of each contribution population
f - float64 fecundity
V_ls - float64 variance at likage equilibrium

Value:
x_prime - float64 density of new juviniles
mu_prime - float64 - mean trait value 
sigma_prime - float64 - variance in trait values
"""
# function reproduction_drift(N, mu, sigma, f, V_le)
#     @assert length(mu) == length(sigma)
#     @assert length(x) == length(mu)
    
#     n = length(mu)
#     epsilon = rand(d,n)
#     mu = mu .+ sigma.*epsilon./sqrt.(N)
    
#     for i in 1:n
#         chi = Distributions.Chisq(N[i]-1)
#         v = rand(chi,1)
#         sigma[i] = sqrt(v*sigma[i]^2/(n-1))
#     end
    
#     mu_prime = sum(mu.*N)/sum(x)
#     G = sum(N .* (sigma.^2 .+ mu.^2))/sum(x) - mu_prime^2
#     sigma_prime = G/2 + V_le/2
#     N_prime = f*sum(N)
#     if sigma_prime < 0
#         print("issue")
#         print("\n")
#         sigma_prime = 0.001
#     end
#     return x_prime, mu_prime, sqrt(sigma_prime) 
# end 

end # module 