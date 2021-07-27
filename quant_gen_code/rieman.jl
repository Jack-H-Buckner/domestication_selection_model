"""
This module contians code to aproximate the distribition of a 
trait within a populaiton and evaluate various operators using 
riemann integrals.

This method is not the most computationally efficent, but it
is reletively straight forward to implement and as such is
a useful way to test the performance of other more efficent
but more complex aproximation methods
"""
module rieman_gen

include("utils.jl")

mutable struct step_function
    a::Float64
    b::Float64
    n::Int
    interval::Float64
    nodes::AbstractVector{Float64}
    values::AbstractVector{Float64}
end

function copy(f)
    return step_function(f.a,f.b,f.n,f.interval,f.nodes,f.values)
end


function init_step_function_unif(a,b,n)
    interval = (b-a)/(n-1)
    nodes = collect(a:interval:b)
    values = repeat([1/interval],n)
    return step_function(a,b,n,interval,nodes,values)
end 


function init_step_function_normal(a,b,n,mu,sigma)
    interval = (b-a)/(n-1)
    nodes = collect(a:interval:b)
    values = broadcast(x -> utils.normal_density(mu, sigma, x),nodes)
    values = values ./ sum(values)
    values = values ./ interval
    return step_function(a,b,n,interval,nodes,values)
end 


function stabalizing_selection!(step_function, x, opt, scale)
    v = step_function.values .* broadcast(x -> utils.squared_exponential(opt,scale,x), step_function.nodes)
    x = x*sum(v)*step_function.interval
    v = v ./sum(v)
    step_function.values = v ./ step_function.interval
    return step_function, x
end 



function stabalizing_selection(step_function1, x, opt, scale)
    v = step_function1.values .* broadcast(x -> utils.squared_exponential(opt,scale,x), step_function1.nodes)
    x = x*sum(v)*step_function1.interval
    v = v ./ sum(v)
    v = v ./ step_function1.interval
    step_function2 = step_function(step_function1.a,step_function1.b,step_function1.n,step_function1.interval,step_function1.nodes,v)
    return step_function2, x
end 


function combine_distributions!(step_function, step_function2, w)
    @assert all(step_function.nodes .== step_function2.nodes)
    v = (1-w) .* step_function.values 
    v = v .+ w .* step_function2.values 
    step_function.values = v
    return step_function
end 



function random_mating!(step_function, Vle)
    
    # update mean 
    difference = (step_function.nodes .+ transpose(step_function.nodes))./2
    product = step_function.values .* transpose(step_function.values)

    v = broadcast(n -> utils.conv(step_function, difference, product, n), 1:step_function.n)
    
    step_function.values = v
    
    # add  Vle to means 
    g = init_step_function_normal(step_function.a,step_function.b,step_function.n,0,sqrt(Vle/2))
    
    difference = g.nodes .+ transpose(step_function.nodes)
    product = g.values .* transpose(step_function.values)
    
    v = broadcast(n->utils.conv_2(step_function,g,product,difference,n), 1:step_function.n)
    
    v = v ./sum(v)
    step_function.values = v ./ step_function.interval
    
    
    return step_function
end 


end # module