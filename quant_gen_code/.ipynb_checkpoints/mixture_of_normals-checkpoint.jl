"""
Many operators used to model the distribution of quantitative traits in a popualtion 
preserve the functional form of normal distibutions and mixtures of normals. This 
is very convienent because quantitative tend to be normally distributed in natural popuatlions.

Two primary operators that are used in these models are stabalizing seleciton and random mating.
stabaliing selection preserves normality and the number of mixture components, and random mating 
preserves normality, but increases the number of mixture components to n(n+1)/2. This is problematic 
because repeated application of the random mating opperator will result in a faster an exponential increase 
in the number of mixture components. This means that despite the existance of a closed form solution
in practice these distributions need to be aproximated.  

One convienent method to aproximate these distibutions is to use mixtures of normals with a fixed 
number of components.  this allows all of the relevant operators to be applied while maintianing a 
colosed form solutoion. The results can then be aproximated using a fixed number of components 
to maintian tractability.
"""
module Mixture_of_Normals
include("utils.jl")
include("normal_trait_distribution.jl")
using NLsolve



"""
This mutable struct stores the data required to evaluate and aproximate
distributions of quantitative traits as mixtures of normal distributions

variances are standard deviaionts!!!!
"""
mutable struct mixture_of_normals
    weights::AbstractVector{Float64} # weights for mixture
    means::AbstractVector{Float64} # centers for each component
    variances::AbstractVector{Float64} # variance for each component, unless aproximation used
    N_components::Int # number of components included in mixture 
    nodes::AbstractVector{Float64}
    sigma::Float64 # variance of mixture components in aproximation 
end 

"""
Evaluates the probability of density of a point x given 
a mixture_of_normals object
"""
function evaluate_density(x, mixture) 
    sum(mixture.weights .* utils.normal_density.(mixture.means, mixture.variances, x))
end

    
"""

"""
function update_nodes!(mixture)
    #mixture.N_components = length(mixture.means)
    min_n = mixture.means[argmin(mixture.means)]
    max_n = mixture.means[argmax(mixture.means)]

    nodes = collect(min_n:((max_n-min_n)/(mixture.N_components-1)):max_n)
    while length(nodes) < mixture.N_components
        nodes = vcat(nodes, [nodes[end] + nodes[end]-nodes[end-1]])
    end 
    mixture.nodes = nodes
    return mixture
end


"""
distance between aproximation and true values
"""
function f_aprox!(fx,x, v, mixture)
    weights = x[2:end]
    sigma = x[1]
    fx[2:end] = v .- utils.normal_density.(transpose(mixture.nodes), sigma, mixture.nodes) * weights
    fx[1] = 1 - sum(weights)
    return fx
end 
    
   
    
"""
jacobian for distance between aproximation and true values
"""
function g_aprox!(gx,x,mixture)
    weights = x[2:end]
    sigma = x[1]
    gx[2:end,2:end] = -1*utils.normal_density.(mixture.nodes, sigma, transpose(mixture.nodes)) 
    gx[1,1] = 0
    gx[1,2:end] .= - 1.0
    gx[2:end,1] = -1 * utils.dsigma_normal_density.(transpose(mixture.nodes),sigma,mixture.nodes) * weights
    return gx
end

"""
takes a mixture of normals object and appoximates the
the distribution with a fixed number of components with
a fixed variance. 
"""
function aproximate_mixture!(mixture)
    mixture = update_nodes!(mixture)
    if mixture.N_components < length(mixture.means)
        n = mixture.N_components
        v = broadcast(x -> evaluate_density(x, mixture), mixture.nodes)    
        f!(fx,x) = f_aprox!(fx,x,v,mixture)
        g!(gx,x) = g_aprox!(gx,x,mixture)
        opt = NLsolve.mcpsolve(f!,g!,zeros(n+1), repeat([Inf], n+1), vcat(mixture.sigma,v), iterations = 1000, ftol = 0.00001)
        #opt = NLsolve.nlsolve(f!,g!,vcat(mixture.sigma,v), iterations = 2000, ftol = 0.001)
        #@assert NLsolve.converged(opt)

        mixture.means = mixture.nodes
        @assert sum(opt.zero[2:end]) < 1.0025
        mixture.weights = opt.zero[2:end]/sum(opt.zero[2:end])
        mixture.variances = repeat([opt.zero[1]], length(mixture.nodes))
    end 
            
    return mixture
end 


    
    
"""
takes a mixture of normals object and appoximates the
the distribution with a fixed number of components with
a fixed variance. 
"""
function aproximate_mixture_sigma!(mixture, sigma)
        
    @assert mixture.N_components <= length(mixture.means)
    v = broadcast(x -> evaluate_density(x, mixture), mixture.nodes)    
    f!(fx,x) = f_aprox!(fx,x,v,mixture)
    g!(gx,x) = g_aprox!(gx,x,mixture)
    opt = NLsolve.nlsolve(f!,g!,vcat(mixture.sigma,v), iterations = 100, tol = 0.0001)
    #@assert NLsolve.converged(opt)

    mixture.means = mixture.nodes
    mixture.weights = opt.zero[2:end]
    mixture.variances = repeat([opt.zero[1]], length(mixture.nodes))
          
    return mixture
end     
    
    
    
    
    
"""
stabalizing selection
applies stabalizing selection oprator to a mixture 
object with trait optimum opt and length scale scale
"""
function stabalizing_selection!(mixture,x,opt,scale)
    for i in 1:length(mixture.weights)
        mixture.weights[i], mixture.means[i],mixture.variances[i]  = 
        normal_trait_distribution.selection_vector(mixture.weights[i], 
        mixture.means[i],mixture.variances[i], opt, scale)
    end
    x = x*sum(mixture.weights)
    mixture.weights = mixture.weights./sum(mixture.weights)
    return mixture, x 
end
    
    
    
"""
stabalizing selection
applies stabalizing selection oprator to a mixture 
object with trait optimum opt and length scale scale
"""
function stabalizing_selection(mixture,x,opt,scale)
    n = length(mixture.weights)
    weights = zeros(n)
    means = zeros(n)
    variances = zeros(n)
    for i in 1:length(mixture.weights)
        weights[i], means[i], variances[i]  = 
        normal_trait_distribution.selection_vector(mixture.weights[i], 
        mixture.means[i],mixture.variances[i], opt, scale)
    end
    x = x*sum(mixture.weights)
    weights = weights./sum(weights)
        
    N = mixture.N_components
    nodes = mixture.nodes
    sigma = mixture.sigma
    return mixture_of_normals(weights,means,variances,N,nodes,sigma), x 
end

function combine_distributions!(mixture,mixture2,w)
    mixture.means = vcat(mixture.means, mixture2.means)
    mixture.variances = vcat(mixture.variances, mixture2.variances) 
    mixture.weights = vcat(mixture.weights,w * mixture2.weights) ./(1+w)
    mixture = update_nodes!(mixture)
    return mixture
end 
    

function random_mating!(mixture,Vle)
    n = length(mixture.means)
    means = (mixture.means .+ transpose(mixture.means)) ./2
    variances = sqrt.((mixture.variances.^2 .+ transpose(mixture.variances).^2 .+ Vle) ./ 2)
    weights = mixture.weights .* transpose(mixture.weights)

    mixture.means = zeros(n*(n+1)÷2)
    mixture.variances = zeros(n*(n+1)÷2)
    mixture.weights = zeros(n*(n+1)÷2)
    k = 0
    for i in 1:n
        for j in 1:n
            if i <= j
                k+=1
                mixture.means[k] = means[i,j]
                mixture.variances[k] = variances[i,j]
 
                if i == j
                    mixture.weights[k] = weights[i,j] 
                else
                    mixture.weights[k] = weights[i,j]  + weights[j,i]
                end
            end
        end
    end
    return mixture
end
    
function mean_var(mixture)
    mean = sum(mixture.means .* mixture.weights)/sum(mixture.weights)
    G = sum(mixture.weights .* (mixture.variances.^2 .+ mixture.means.^2))/sum(mixture.weights) - mean^2
    return mean, sqrt(G)
end 
    
end # module