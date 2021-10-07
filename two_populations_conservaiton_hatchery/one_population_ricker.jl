module one_populaiton_ricker

mutable struct population
    N::Float64
    z::Float64
    r::Float64
    K::Float64
    s::Float64
    z1::Float64
end 


function update_populaiton!(model::population,m::Float64)
    N = model.N
    r = model.r
    K = model.K
    s = model.s
    z1 = model.z1
    
    model.N = N*exp(r*(1-N/K) - 0.5*s/(1+s)*model.z^2)
    model.z = N/(N+m)* model.z/(1+s) + z1*m/(N+m)
    model.N += m
    
    return model
end


function equilibrium_population_fpi(model::population,m::Float64)
    tol = 0.000001
    
    state_0 = vcat(model.N,model.z)
    model = update_populaiton!(model,m)
    state_1 = vcat(model.N,model.z)
    i = 0
    while any(abs.(state_0 .- state_1) .> repeat([tol],2))&& i < 100000
        
        i += 1
        state_0 = state_1
        model = update_populaiton!(model,m)
        state_1 = vcat(model.N,model.z)
        
    end
    return model 
end 



end # module 