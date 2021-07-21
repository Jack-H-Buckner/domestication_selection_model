module two_populations_conservaiton_hatchery
include("normal_trait_distribution.jl")
"""
Mutatable struct that stores states
and parameters for the model. 
"""
mutable struct population
    # states 
    x::AbstractVector{Float64} # length 2 - populaiton sizes 
    mu::AbstractVector{Float64} # length 2 - mean trait values
    sigma::AbstractVector{Float64} # length 2 - trait variances 
    # parameters 
    mu_s::AbstractVector{Float64} # length 3 - optimal trait values 
    sigma_s::AbstractVector{Float64} # length 3 - strength of selection
    k::AbstractVector{Float64} # length 2 - density dependence 
    f::AbstractVector{Float64} # length 2 - fecundity 
    V_le::Float64
    s::Float64 # dispersal 
    prop::Float64 # brood stock
    max::Float64 # brood stock 
end


"""
update_populations(x,mu,sigma,r)
updated popautlion dynamics after one time step

Arguments:
x - vector length 2 size of population right after reproduction
mu - vector length 2 mean trait values after reporduction
sigma - vector length 2 variacne of trait values after reproduction 
r - float64 number of individuals released 

Value:

"""
function update_populations(model,releases)
    # add releases to vector 
    x = [model.x[1],releases,model.x[2]]
    mu = [model.mu[1],model.mu[1],model.mu[2]]
    sigma = [model.sigma[1],model.sigma[1],model.sigma[2]]
    
    x[2] = releases
    
    # selection 
    x[1], mu[1], sigma[1] = normal_trait_distribution.selection(x[1], mu[1], sigma[1], model.mu_s[1], model.sigma_s[1])
    v, mu[2], sigma[2] = normal_trait_distribution.selection(x[2], mu[2], sigma[2], model.mu_s[2], model.sigma_s[2])
    x[3], mu[3], sigma[3] = normal_trait_distribution.selection(x[3], mu[3], sigma[3], model.mu_s[3], model.sigma_s[3])
    
    # dispersal
    x1 = [x[1],(1-model.s)*x[2]] 
    x2 = [x[3], model.s*x[2]]
 
    # density dependence
    x1 = normal_trait_distribution.competition(x1,model.f[1],model.k[1])
    x2 = normal_trait_distribution.competition(x2,model.f[2],model.k[2])
    
    # get brood stock
    x1[1], b = normal_trait_distribution.brood_stock(x1[1],model.prop,model.max)
   
    
    # reproduction 
    x1, mu1, sigma1 = normal_trait_distribution.reproduction(x1, [mu[1],mu[2]], [sigma[1],sigma[2]], model.f[1], model.V_le)
    x2, mu2, sigma2 = normal_trait_distribution.reproduction(x2, [mu[3],mu[2]], [sigma[3],sigma[2]], model.f[2], model.V_le)
    
    # update model object
    model.x = [x1,x2]
    model.mu = [mu1,mu2]
    model.sigma = [sigma1,sigma2]
    
    return model
end



function equilibrium_population_fpi(model,releases)
    tol = 0.0001
    
    state_0 = vcat(model.x,model.mu,model.sigma)
    model = update_populations(model,releases)
    state_1 = vcat(model.x,model.mu,model.sigma)
    i = 0
    while any(abs.(state_0 .- state_1) .> repeat([tol],6))&& i < 100000
        
        i += 1
        state_0 = state_1
        model = update_populations(model,releases)
        state_1 = vcat(model.x,model.mu,model.sigma)
        
    end
    return model 
end 





function update_populations_return_state(model,releases)
    # add releases to vector 
    x = [model.x[1],releases,model.x[2]]
    mu = [model.mu[1],model.mu[1],model.mu[2]]
    sigma = [model.sigma[1],model.sigma[1],model.sigma[2]]
    
    x[2] = releases
    
    # selection 
    x[1], mu[1], sigma[1] = normal_trait_distribution.selection(x[1], mu[1], sigma[1], model.mu_s[1], model.sigma_s[1])
    v, mu[2], sigma[2] = normal_trait_distribution.selection(x[2], mu[2], sigma[2], model.mu_s[2], model.sigma_s[2])
    x[3], mu[3], sigma[3] = normal_trait_distribution.selection(x[3], mu[3], sigma[3], model.mu_s[3], model.sigma_s[3])
    
    x_out = x
    mu_out = mu
    sigma_out = sigma
    # dispersal
    x1 = [x[1],(1-model.s)*x[2]] 
    x2 = [x[3], model.s*x[2]]
     
    
    # density dependence
    x1 = normal_trait_distribution.competition(x1,model.f[1],model.k[1])
    x2 = normal_trait_distribution.competition(x2,model.f[2],model.k[2])
    
    # get brood stock
    x1[1], b = normal_trait_distribution.brood_stock(x1[1],model.prop,model.max)
   
    
    # reproduction 
    x1, mu1, sigma1 = normal_trait_distribution.reproduction(x1, [mu[1],mu[2]], [sigma[1],sigma[2]], model.f[1], model.V_le)
    x2, mu2, sigma2 = normal_trait_distribution.reproduction(x2, [mu[3],mu[2]], [sigma[3],sigma[2]], model.f[2], model.V_le)
    
    # update model object
    model.x = [x1,x2]
    model.mu = [mu1,mu2]
    model.sigma = [sigma1,sigma2]
    
    return model, x_out, mu_out, sigma_out
end







end # module