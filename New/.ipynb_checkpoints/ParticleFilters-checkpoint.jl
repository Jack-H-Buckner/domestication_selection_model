"""
some function to use the particle filtering algorithm to compute likelihood functions
for model parameters. 
"""
module ParticleFilters 

using StatsBase
using LogExpFunctions
using Distributions
include("pedigreeAproximation.jl")
include("AgeTraitStructuredModels.jl")

#######################
### Basic machinery ###
#######################

mutable struct particleFilter
    N::Int64
    logLik::AbstractVector{Float64}
    particles#::AbstractVector{pedigreeAproximation.population}
end


function overwrite_state!(pop_replace,pop_copy)
    pop_replace.abundanceN .= pop_copy.abundanceN
    pop_replace.abundanceH .= pop_copy.abundanceH
    pop_replace.zH .= pop_copy.zH
    pop_replace.pi .= pop_copy.pi
    pop_replace.z .= pop_copy.z
end 



"""
Resamples popualtion objects in the particle filter according to the weights. Assumes
the parameter values are equal across all population objects and only states need to
be updated. 
"""
function resample_particle_state!(particleFilter)
    
    # sample indecies to save
    inds = sample(1:particleFilter.N, weights(softmax(particleFilter.logLik)),particleFilter.N, replace = true)
    # create list of indecies to save
    unique_inds = unique(inds)
    # creat list of number of times to copy each popualtion
    N_inds = zeros(length(unique_inds))
        
    for i in 1:length(unique_inds) N_inds[i] = sum(unique_inds[i] .== inds) end 
    
    # create a list of indecies that can be over writen 
    overwrite_inds = collect(1:particleFilter.N)[broadcast(ind -> !(ind in unique_inds), 1:particleFilter.N)]
    
    # over write 
    i = 0 # location in unique_inds, N_inds
    j = 0 # location in overwrite_inds
    for ind_copy in unique_inds
        i += 1
        for k in 2:N_inds[i]
            j += 1
            ind_overwrite = overwrite_inds[j] 
            overwrite_state!(particleFilter.particles[ind_overwrite],particleFilter.particles[ind_copy])
        end 
    end 
    particleFilter.logLik = zeros(particleFilter.N)
end



############################
### Likelihood functions ###
############################



function recruits_pedigree(pedigrees,population)
    min_p = 10^(-30)
    # map pedigrees to grid
    pi_ls = zeros(length(pedigrees))
    i = 0  
    for p in pedigrees
        i += 1
        pi_ls[i] = round.(p .* 2^population.n)
            
    end 
        
    lL = 0.0

    for p in pi_ls
        lL += log(population.pi[floor(Int,p)+1,1] + min_p)
    end 
    return lL
end 


"""
log noraml o
"""
function recruitsN(estimate,precision,population)
    N = log(population.abundanceN[1] + 10^(-30))
        
    return log(sqrt(precision)/sqrt(2*pi*estimate^2)) - 0.5*precision*(log(estimate) -N)^2
end 

"""
Binomial likelihood for natrual/ hatchery recruits 
"""
function recruitsH(estimate,precision,population)
    N = log(population.abundanceH[1] + 10^(-30))
        
    return log(sqrt(precision)/sqrt(2*pi*estimate^2)) - 0.5*precision*(log(estimate) -N)^2
end 


#############################
## Model spcific functions ##
#############################

#### parameters ####

# T_start - time when stocking begins - (known)
# agestructure - demographic parameters/ rates  - (known)
# Nt releases - number of hathceyr fish releaseed - (known)

# p_im - probability of an indiviudal immigating - (estimated)
# s - strength of selection reletive to genetic variance - (estimated)
# RRS - reletive reproductive success/ fitness - (estimated)
# p_spawn - number of migrants that spawn - (optional)

#### data ####
# t_pedigree  - time steps of pedigree estiamtes  
# t_recruitsN - time steps of recruitment estiamtes 
# t_recruitsH
# pedigrees - pedigree data 
# recruitsN  - recruits data 
# recruitsH
# precision - precision of recruitment estiamtes 

#### MCMC parameters ####
# N_particles


struct data
    T_end::Int64 # time point of last data point 
    ageStructure # age stucture parameters 
    Nt::AbstractVector{Int} # number of hatchery releases each year 
    time_step!::Function 
    t_pedigree::AbstractVector{Int} 
    t_recruitsN::AbstractVector{Int}
    t_recruitsH::AbstractVector{Int} 
    pedigrees::AbstractVector{AbstractVector{Float64}}
    recruitsN::AbstractVector{Float64}
    recruitsH::AbstractVector{Float64}
    precision::Float64
end 



mutable struct params
    p_im::Float64
    s::Float64
    RRS::Float64
    p_spawn::Float64
end 

# convert RRS to selection strength and trait difference 

function likelihood(params, data; N_particles = 100, n = 5)
    # Convert RRS to mu_im
    mu_im = AgeTraitStructuredModels.solve_trait_difference(params.RRS,params.s)

    # initialize particles
    logLik = 0.0
    logLikls = zeros(N_particles)
    populaitons = []
    for i in 1:N_particles
        push!(populaitons, pedigreeAproximation.init_population(data.ageStructure,  0.0, params.s, n, params.p_spawn))
            
    end 
    pf = particleFilter(N_particles,logLikls,populaitons)
    
    t_data = unique(vcat(data.t_pedigree,data.t_recruitsN,data.t_recruitsH))
        
    for t in 1:data.T_end
   
        # update particles
        # 
        Threads.@threads for i in 1:N_particles
            
            if data.Nt[t] > 0
                data.time_step!(populaitons[i], rand(Distributions.Binomial(data.Nt[t], params.p_im)), mu_im) 
            else
                data.time_step!(populaitons[i], 0, mu_im) 
            end
        

            ## pedigress
            if t in data.t_pedigree
                pf.logLik[i] += recruits_pedigree(data.pedigrees[data.t_pedigree .== t][1],pf.particles[i])
            end 


            # natrual recruitment
            if t in data.t_recruitsN
                pf.logLik[i] += recruitsH(data.recruitsN[data.t_recruitsN .== t][1],data.precision,pf.particles[i])
            end 



            # hatchery migation
            if t in data.t_recruitsH
                pf.logLik[i] += recruitsH(data.recruitsH[data.t_recruitsH .== t][1],data.precision,pf.particles[i])
            end 
  
        end 
        
        # likelihood estiamte using log sum exponential trick
        maxLL = pf.logLik[argmax(pf.logLik)]
        scaled_weights = exp.(pf.logLik .- maxLL)
        mean_scaled_weight = sum(scaled_weights/pf.N) 
        logLik += log(mean_scaled_weight*exp(maxLL))

        # resample particles 
        resample_particle_state!(pf) 
        
    end 
    
    return logLik
    
end 



######################################
### Metropolis Hastings algorithm ###
######################################

"""
    A functor that stores data to define a proposal distirbution 
on a bounded interval. 
"""
struct sampler
    dims::Int64
    a::AbstractVector{Float64}
    b::AbstractVector{Float64}
    c::AbstractVector{Float64}
end 
"""
    returns a new sample yt drawn form a beta distribution in each dimension
and a correction factor cor to account for a asymetry of the distribution. 
"""
function (p::sampler)(xt)
    # rescale xt to unit interval
    xt = (xt .- p.a) ./ (p.b.-p.a) 
    yt = zeros(p.dims)
    cor = 1.0
    for i in 1:p.dims
        # define p(Yt|xt) sample yt and evaluate p(yt|xt)
        alpha = xt[i]*p.c[i]
        beta = (1-xt[i])*p.c[i]
        d = Distributions.Beta(alpha + 1.0,beta+1.0)
        yt[i] = rand(d)
        qyx = pdf(d,yt[i])
        
        # evaluate p(xt|yt)
        alpha = yt[i]*p.c[i]
        beta = (1-yt[i])*p.c[i]
        d = Distributions.Beta(alpha+1.0,beta+1.0)
        qxy = pdf(d,xt[i])
        
        # compute correction (in this dimension)
        cor *= qxy/qyx
    end
    return (yt .* (p.b.-p.a)) .+ p.a , cor
end 

function MCMC(x0, sampler, data; Nmc = 2000)
    xt, cor = sampler(x0)
    pars = params(xt[1],xt[2],xt[3], 1.0)
    samples = zeros(sampler.dims, Nmc)
    lLx = likelihood(pars, data;N_particles = 10, n = 4)
    
    for i in 1:Nmc
        if mod(i, 10) == 0 print(i, " ") end 
        yt, cor = sampler(xt)
        pars = params(yt[1],yt[2],yt[3], 1.0)
        lLy = likelihood(pars, data)
        pars = params(xt[1],xt[2],xt[3], 1.0)
        lLx = likelihood(pars, data)
        p = cor*exp(lLy - lLx)
        if p > rand()
            xt = yt
            lLx = lLy
        end
        samples[:,i] = xt
    end 
    
    return samples
end


end 