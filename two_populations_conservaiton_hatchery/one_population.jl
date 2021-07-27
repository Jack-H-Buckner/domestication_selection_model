module one_population_model
include("../quant_gen_code/normal_trait_distribution.jl")
include("../quant_gen_code/mixture_of_normals.jl")
include("../quant_gen_code/rieman.jl")

## normals aproximation 

mutable struct population 
    x::Float64
    mu::Float64
    sigma::Float64
    opt::Float64
    scale::Float64
    Vle::Float64
    r::Float64
    k::Float64
end



function update_populaion!(model)
    x = model.x
    mu = model.mu
    sigma = model.sigma
    r = model.r
    k = model.k
    Vle = model.Vle
    opt = model.opt
    scale = model.scale
    
    x, mu, sigma = normal_trait_distribution.reproduction(x, mu, sigma, r, Vle)
    
    x, mu, sigma = normal_trait_distribution.selection(x, mu, sigma, opt, scale)
    
    x = normal_trait_distribution.competition(x,r,k)
    model.x = x
    model.mu = mu
    model.sigma = sigma
    return model 
end 




mutable struct population_hatchery
    x::Float64
    mu::Float64
    sigma::Float64
    opt::Float64
    scale::Float64
    opt_H::Float64
    scale_H::Float64
    Vle::Float64
    r::Float64
    k::Float64
    prop::Float64
    max::Float64
end



function update_populaion_hatchery!(model, releases)
    x = model.x
    mu1 = model.mu
    sigma1 = model.sigma
    r = model.r
    k = model.k
    Vle = model.Vle
    opt = model.opt
    scale = model.scale
    opt_H = model.opt_H
    scale_H = model.scale_H
    
    # selection

    x, mu, sigma = normal_trait_distribution.selection(x, mu1, sigma1, opt, scale)
    v, mu_H, sigma_H = normal_trait_distribution.selection(x, mu1, sigma1, opt_H, scale_H)

    
    
    x = [x,releases]
    
    # densiy dependence
    x = normal_trait_distribution.competition(x,r,k)
    
    
    
    # brood stock
    x[1], b = normal_trait_distribution.brood_stock(x[1],model.prop,model.max)

    # reproduction 
    x, mu, sigma = normal_trait_distribution.reproduction(x, [mu,mu_H], [sigma,sigma_H], r, Vle)
    

    
    model.x = x
    model.mu = mu
    model.sigma = sigma
    
    return model 
end


function equilibrium_population_fpi(model,releases)
    tol = 0.000001
    
    state_0 = vcat(model.x,model.mu,model.sigma)
    model = update_populaion_hatchery!(model,releases)
    state_1 = vcat(model.x,model.mu,model.sigma)
    i = 0
    while any(abs.(state_0 .- state_1) .> repeat([tol],3))&& i < 100000
        
        i += 1
        state_0 = state_1
        model = update_populaion_hatchery!(model,releases)
        state_1 = vcat(model.x,model.mu,model.sigma)
        
    end
    return model 
end 



### mixture or normals aproximation


mutable struct population_hatchery_mixture
    x::Float64
    dsn
    opt::Float64
    scale::Float64
    opt_H::Float64
    scale_H::Float64
    Vle::Float64
    r::Float64
    k::Float64
    prop::Float64
    max::Float64
end



function update_populaion_hatchery_mixture!(model, releases)
    x = model.x
    r = model.r
    k = model.k
    Vle = model.Vle
    opt = model.opt
    scale = model.scale
    opt_H = model.opt_H
    scale_H = model.scale_H
    
    dsn = model.dsn

    
    dsn, x = Mixture_of_Normals.stabalizing_selection!(dsn,x,opt,scale)
    dsn_h, x_h = Mixture_of_Normals.stabalizing_selection(dsn,x,opt_H,scale_H)
    

    # combine distribuitons 
    dsn = Mixture_of_Normals.combine_distributions!(dsn,dsn_h,releases/x)

    x = [x,releases]
    x = normal_trait_distribution.competition(x,r,k)
    
    

    dsn = Mixture_of_Normals.random_mating!(dsn,Vle)


    dsn = Mixture_of_Normals.aproximate_mixture!(dsn)

    model.dsn = dsn
    model.x = r*sum(x)
    return model
end


function equilibrium_population_fpi_mixture!(model,releases)
    tol = 0.0001

    state_0 = model.x
    model = update_populaion_hatchery_mixture!(model,releases)
    state_1 = model.x
    i = 0
    while (abs(state_0 .- state_1) > tol)&& i < 100000

        i += 1
        state_0 = state_1
        model = update_populaion_hatchery_mixture!(model,releases)
        state_1 = model.x
        
    end
    return model 
end 





### step function + rieman integral aproximation 


mutable struct population_hatchery_rieman
    x::Float64
    dsn
    opt::Float64
    scale::Float64
    opt_H::Float64
    scale_H::Float64
    Vle::Float64
    r::Float64
    k::Float64
    prop::Float64
    max::Float64
end



function update_populaion_hatchery_rieman!(model, releases)
    x = model.x
    r = model.r
    k = model.k
    Vle = model.Vle
    opt = model.opt
    scale = model.scale
    opt_H = model.opt_H
    scale_H = model.scale_H
    
    
    dsn1 = rieman_gen.copy(model.dsn)

    
    dsn_h, x_h = rieman_gen.stabalizing_selection(dsn1,x,opt_H,scale_H)
    dsn, x = rieman_gen.stabalizing_selection!(model.dsn,x,opt,scale)
    

  
    # combine distribuitons 
    
    dsn = rieman_gen.combine_distributions!(dsn,dsn_h, releases/(releases + x))
    
    
    
    x = [x,releases]
    x = normal_trait_distribution.competition(x,r,k)
    
    
    
    # brood stock
    x[1], b = normal_trait_distribution.brood_stock(x[1],model.prop,model.max)
    
 

    #reproduction
    x = r*sum(x)
    dsn = rieman_gen.random_mating!(model.dsn,Vle)
    

    model.dsn = dsn
    model.x = sum(x)
    
    return model
end


function equilibrium_population_rieman_fpi(model,releases)
    tol = 0.00001
    
    state_0 = model.x
    model = update_populaion_hatchery_rieman!(model,releases)
    state_1 = model.x
    i = 0
    while any(abs.(state_0 .- state_1) .> repeat([tol],1))&& i < 100000
        
        i += 1
        state_0 = state_1
        model = update_populaion_hatchery_rieman!(model,releases)
        state_1 = model.x
        
    end
    return model 
end 



end # module