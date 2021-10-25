"""
This module contains code for a simple individual based model that 
describes the effects of gene flow and seleciton on the distribution of 
a trait in a populaiton. 

It also tracks the proportion of each individuals ancestry in the original population 
and that immigated to the population. The model also trackes the population of origin of
each indiviudals parents and grand parents. 

The goal of this code is to develop some intuition for the information on the
fittness effects of gene-flow contained in data on pedigrees. 

The module contains two data structures to define individuals and populaitons 
and a set of function to update the population. 
"""
module IBM_pedigrees

"""
Data structure for individuals in the populaiton
"""
struct individual
    g::Float64
    p::Float64
    t_macro
    t_refined
end 

"""
Data structure for the population
"""
mutable struct population
    d
    list::AbstractVector{individual}
end 

"""
Initializes a populaiton with N individuals 
Vle is the variance of trait distribution at equibrium 
"""
function init_population(N,Vle)
    d = Distributions.Normal(0,Vle/2)
    g = rand(d,N)
    p = zeros(N)
    list = []
    for i in 1:N
        list = vcat(list,individual(g[i],p[i], "W", "WW"))
    end 
    return population(d,list)
end 

    

function micro_to_macro(state)
    if state == "H1H1"
        return "H" 
    elseif occursin("H1", state)
        return "F1"
    else
        return "W"
    end 
end 

function reduce_to_micro(state)
    if state == "WH1"
        return "H1W" 
    elseif state == "WF1"
        return "F1W"
    elseif state == "HF1"
        return "F1H"
    elseif state == "HW"
        return "WH"
    elseif state == "HH1"
        return "H1H"
    elseif state == "F1H1"
        return "H1F1"
    else
        return state
    end 
end 

    

function reproduce(population)
    pair = rand(population.list,2)
    g = 0.5*pair[1].g + 0.5*pair[2].g + rand(population.d,1)[1]
    p = 0.5*pair[1].p+0.5*pair[2].p
    t_micro = reduce_to_micro(pair[1].t_macro * pair[2].t_macro)
    t_macro = micro_to_macro(t_micro)
    ind = individual(g,p,t_macro,t_micro)
    return ind
end 


function reproduce!(population,N)
    list = []
    for i in 1:N
        list = vcat(list,reproduce(population))
    end 
    population.list = list
    return population
end 

function selection!(population,s)
    p= zeros(length(population.list))
    w= zeros(length(population.list))
    i = 0
    for ind in population.list
        i += 1
        p[i] = rand()[1]
        w[i] = exp(-s/2*ind.g^2)
    end 
    population.list = population.list[p.<w]
    return population
end

function immigration!(population,M,g_im)
    g = rand(population.d,M) .+ g_im
    p = ones(M)

    list = []
    for i in 1:M
        list = vcat(list,individual(g[i],p[i], "H1", "H1"))
    end 
    population.list = vcat(population.list,list)
    return population
end
        
function distributions(population)
            
    N = length(population.list)
    g = zeros(N)
    p = zeros(N)

    for i in 1:N
        ind = population.list[i]
        g[i] = ind.g
        p[i] = ind.p
    end 
    return g,p
end
        
function data(population, path)
    N = length(population.list)
    g = zeros(N)
    p = zeros(N)
    t_micro = []
    t_macro = []
    for i in 1:N
        ind = population.list[i]
        g[i] = ind.g
        p[i] = ind.p
        t_micro = vcat(t_micro,ind.t_refined)
        t_macro = vcat(t_macro,"a"*ind.t_macro)
    end 
    df = DataFrames.DataFrame(g_=g, p_=p, 
            micro_=t_micro, macro_=t_macro)
   CSV.write(path, df)
end 

end # module 