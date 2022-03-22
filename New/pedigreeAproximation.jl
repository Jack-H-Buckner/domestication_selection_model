module pedigreeAproximation
    
include("AgeStructuredModels.jl")
using Distributions

mutable struct population
    # aproximation meta data
    n::Int64 # grid refinement 
    m::Int64 # length of vectors 2^n + 1
    
    # genetic state
    pi::AbstractMatrix{Float64} # pedigree dsn natrual 
    z::AbstractMatrix{Float64} # trait values natural 
    V::Float64 # genetic variance conditional on pedigree 
    
    # demgraphic state
    abundance::AbstractVector{Float64}
    
    # parameters
    ageStructure 
    correction::Float64
    s::Float64 # selection strength 
    theta::Float64 # optimal trait value 
    
end 


mutable struct immigrants
    # aproximation meta-data
    n::Int64
    m::Int64
    mean::Float64
    V::Float64
    
    # genetic state
    N::Float64
    pi::AbstractVector{Float64} # pedigree dsn natrual 
    z::AbstractVector{Float64} # trait values natural 
end 



"""
     V_star_prime(Vle, s)

Computes equilibrium variance given severgation varinace Vle/2 and
slection strength s. measurement after reproduciton but before selection
"""
function V_star_prime(Vle, s)
    #sigma_s = 1/s
    V_prime = V -> (1/V + s)^(-1)/2 + Vle/2#V -> (V*sigma_s^2/(V + sigma_s^2))/2 + Vle/2
    V0 = Vle
    V1 = V_prime(Vle)
    while abs(V0-V1) > 10^-6
        V0 = V1
        V1 = V_prime(V0)
    end 
    return V1
end 


function correction(Vle, s, theta)
    
    # set trait at optimum value and variance = Vle
    grid = collect(-30:0.01:30)
    d = Distributions.Normal(theta, sqrt(V_star_prime(Vle, s)))
    trait = pdf.(d, grid)
    trait = trait ./ sum(trait)
    correction = 1/sum(trait .* exp.(-s/2 .* (grid .- theta).^2))
    
    return correction
end


function init_population(ageStructure,  theta, s,n)
    # aproximation meta data
    m = 2^n + 1
    
    # genetic state
    pi = zeros(m);pi[1] = 1.0
    z = zeros(m); z[1] = theta
    V = V_star_prime(1.0, s)
    
    pi = transpose(repeat(transpose(pi),ageStructure.Amax))
        
    z = transpose(repeat(transpose(z),ageStructure.Amax))
        
    # parameters
    correct = correction(1.0, s, theta)
    
    # demographic state
    abundance = AgeStructuredModels.stable_age_structure(ageStructure)
    
    return population(floor(Int,n), floor(Int, m),pi,z,V,abundance,ageStructure,correct,s,theta)
end 


function init_immigrants(population,N,mean)
    n = population.n
    m = population.m
    pi = zeros(m); pi[end] = 1.0
    z = zeros(m); z[end] = mean

    return immigrants(floor(Int,n),floor(Int,m),mean,population.V,N,pi,z)
end 


function pedigree_trait_update!(z, pi)
    pi_mat = pi .* transpose(pi)
    pi = round_conv!(pi,pi_mat)
    z_mat  = (z .+ transpose(z))./2
    z_mat =  pi_mat .* z_mat 
    round_conv!(z,z_mat)
    z ./= pi # .+ 10^(-50))
    for i in 1:length(z)
        if isnan(z[i])
            z[i] = 0
        end 
    end 
    return z, pi
end


function round_conv!(v,mat)
    # calcualte length for loops 
    n = size(mat)[1]
    
    # zero out v
    for i in 1:n v[i] = 0 end 
    
    #compute concolution
    for i in 1:n
        for j in 1:n
            if mod(i+j,2) == 0
                v[floor(Int,(i+j)/2)] += mat[i,j]
            else
                v[floor(Int,(i+j)/2)] += 0.5*mat[i,j]
                v[floor(Int,(i+j)/2)+1] += 0.5*mat[i,j]
            end 

        end
    end
    return v
end 


function reproduction(population)

    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    
    if f_total == 0
        zeros(population.m),zeros(population.pi),0.0
    end 

    pi = population.pi * (population.abundance .*  population.ageStructure.Fecundity) ./ f_total

    
    z = population.pi .* population.z * (population.abundance .*  population.ageStructure.Fecundity) 
    

    z = z ./ (population.pi * (population.abundance .*  population.ageStructure.Fecundity) .+ 10^(-30))
    
    # convolution - random mating 
    z, pi = pedigree_trait_update!(z, pi)
     
    return z, pi, f_total
end 

function prop_survival(z,z_prime,s,theta,V, V_prime)
    p1 = sqrt(V_prime)/sqrt(V ) 
    p2 = exp(-1/2 *( (s*theta^2 + z^2/V) - z_prime^2/V_prime))
    return p1*p2
end 


function sampleN(n,pi)
    if pi < 10^(-6)
        return 0
    elseif n < 10^(-6)
        return 0
    else
        return rand(Distributions.Binomial(n, pi))
            
    end
end 

function selection(z,pi,N,population)
    theta = population.theta
    s = population.s
    V = population.V
    V_prime = (s+1/V)^(-1)
    z_prime = V_prime.*(s*theta .+ z./V)
    p_survival = broadcast(i -> prop_survival(z[i],z_prime[i],s,theta,V, V_prime),1:length(z))
 
    Nvec = broadcast(n -> floor(Int, n), N.*pi)
    Nvec = broadcast(i->sampleN(Nvec[i],p_survival[i]), 1:population.m)
    if sum(Nvec) == 0
        return zeros(population.m), zeros(population.m), sum(Nvec)
    end 
    pi = Nvec ./ sum(Nvec)
    return z_prime, pi, sum(Nvec)
end 


function immigration(z,pi,N, immigrants)
    Nnew = N + immigrants.N
    
    # allow for demographic stochasticity in trait value around mean 
    immigrants.z[end] = immigrants.mean
    
    # new proportions and means are averages
    pi = (N*pi .+ immigrants.N*immigrants.pi)./Nnew
    
    z = (N*pi.*z + immigrants.N*immigrants.pi.*immigrants.z ) ./ (N*pi .+ immigrants.N*immigrants.pi .+ 10^(-30))
    
    
    return z, pi, Nnew
end 

function sampleR(y)

    if y < 10^(-6)
        return 0
    else
        
        return rand(Distributions.Poisson(y))
            
    end 
end 



"""
Beverton Holt recruitment 

f_total - number of eggs
populations - populaiton objct for parameters
before - boolian does selection occur before or after selection 
"""
function recruitment(z, pi, N, population, before)
    
    if before 
        N *= population.correction 
        R = population.ageStructure.SRCurve(N)
    else
        R = population.correction * population.ageStructure.SRCurve(N)
    end 

    y = R*pi
    Nvec = broadcast(x->sampleR(x),y)
            
    pi = Nvec./sum(Nvec)
    return z, pi, sum(Nvec)
end 


"""
Updates the age structure of the popuatlion and adds recruits
"""
function ageing!(population, z, pi, N)

    
    new_z = zeros(population.m,population.ageStructure.Amax)
    new_z[:,2:end] = population.z[:,1:end-1]
    new_z[:,1] = z
    
    new_pi = zeros(population.m,population.ageStructure.Amax)
    new_N = zeros(population.ageStructure.Amax)
    for a in 2:population.ageStructure.Amax
        
        
        pi_a = population.pi[:,a-1]
        Npi_a = population.abundance[a-1] .* pi_a
        sa = population.ageStructure.Survival[a-1]
       
        Npi_a = broadcast(n -> sampleN(round(Int,n),sa),Npi_a)
        
        if sum(Npi_a) > 0
            new_pi[:,a] = Npi_a / sum(Npi_a)
        else
            new_pi[:,a] = zeros(population.m)
                
        end 
        #new_pi[:,a] = pi_a 
        new_N[a] = sum(Npi_a)#population.abundance[a-1] * sa
    end 
    
    new_pi[:,1] = pi
    new_N[1] = N 
 
    population.abundance = new_N
    population.z= new_z
    population.pi= new_pi
    
end 



function time_step_DSI!(population)
    z, pi, N = reproduction(population)
    z, pi, N =  recruitment(z, pi, N, population, true)
    z, pi, N = selection(z,pi,N,population)
    ageing!(population, z, pi, N)
end 


function time_step_DSI!(population, immigrants)
    z, pi, N = reproduction(population)
    z1, pi1, N1 =  recruitment(z, pi, N, population, true)
    z1, pi1, N1 = selection(z1,pi1,N1,population)
    z, pi, N = immigration(z1,pi1,N1, immigrants)
    ageing!(population, z, pi, N)
end




function SSB(population)
    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    return f_total ./ population.ageStructure.Fecundity[argmax(population.ageStructure.Fecundity)]
end 
function Recruits(populaiton)
    return populaiton.abundance[1]
end 

function trait_moments(population)
    pi = population.pi[:,1]
    z = population.z[:,1]
    return sum(pi.*z) 
end 



function fitness(population)
    
    z,pi,N = reproduction(population)
    theta = population.theta
    s = population.s
    V = population.V
    V_prime = (s+1/V)^(-1)
    z_prime = V_prime.*(s*theta .+ z./V)
    p_survival = broadcast(i -> prop_survival(z[i],z_prime[i],s,theta,V, V_prime),1:length(z))
    

    return sum(p_survival .* pi)
    
end 


function distribution(x,population)
    pi = population.pi[:,1]
    z = population.z[:,1]
    v = population.V
    return sum(broadcast(i -> pi[i]*pdf(Distributions.Normal(z[i],sqrt(v)),x), 1:length(pi)))
end 





##################################################
### define new methods that only track pedigrees #
##################################################


mutable struct population_fixed
    # aproximation meta data
    n::Int64 # grid refinement 
    m::Int64 # length of vectors 2^n + 1
    
    # genetic state
    pi::AbstractMatrix{Float64} # pedigree dsn natrual 
    z::AbstractVector{Float64} # trait values natural 
    V::Float64 # genetic variance conditional on pedigree 
    
    # demgraphic state
    abundance::AbstractVector{Float64}
    
    # parameters
    ageStructure 
    correction::Float64
    s::Float64 # selection strength 
    theta::Float64 # optimal trait value 
    
end 

function init_population_fixed(ageStructure, theta, s,n, mean)
    # aproximation meta data
    m = 2^n + 1
    
    # genetic state
    pi = zeros(m);pi[1] = 1.0
    z = collect(theta:((mean -theta)/2^n):mean)

    V = V_star_prime(1.0, s)
    
    pi = transpose(repeat(transpose(pi),ageStructure.Amax))
        
    # parameters
    correct = correction(1.0, s, theta)
    
    # demographic state
    abundance = AgeStructuredModels.stable_age_structure(ageStructure)
    
    return population_fixed(floor(Int,n), floor(Int, m),pi,z,V,abundance,ageStructure,correct,s,theta)
end 

function reproduction_fixed(population)

    f_total = sum(population.abundance .* population.ageStructure.Fecundity)
    
    if f_total == 0
        zeros(population.m),zeros(population.pi),0.0
    end 
    pi = population.pi * (population.abundance .*  population.ageStructure.Fecundity) ./ f_total
    pi_mat = pi .* transpose(pi)
    pi = round_conv!(pi,pi_mat)
     
    return pi, f_total
end 



function selection_fixed(pi,N,population)
    theta = population.theta
    s = population.s
    V = population.V
    z = population.z
    V_prime = (s+1/V)^(-1)
    z_prime = V_prime.*(s*theta .+ z./V)
    p_survival = broadcast(i -> prop_survival(z[i],z_prime[i],s,theta,V, V_prime),1:length(z))
 
    Nvec = broadcast(n -> floor(Int, n), N.*pi)
    Nvec = broadcast(i->sampleN(Nvec[i],p_survival[i]), 1:population.m)
    if sum(Nvec) == 0
        return zeros(population.m), zeros(population.m), sum(Nvec)
    end 
    pi = Nvec ./ sum(Nvec)
    return pi, sum(Nvec)
end 



function recruitment_fixed(pi, N, population, before)
    
    if before 
        N *= population.correction 
        R = population.ageStructure.SRCurve(N)
    else
        R = population.correction * population.ageStructure.SRCurve(N)
    end 

    y = R*pi
    Nvec = broadcast(x->sampleR(x),y)
            
    pi = Nvec./sum(Nvec)
    return pi, sum(Nvec)
end 



function immigration_fixed(pi,N, N_im)
    Nnew = N + N_im
    pi = N.*pi./Nnew
    pi[end] += N_im / Nnew
    return pi, Nnew
end 




"""
Updates the age structure of the popuatlion and adds recruits
"""
function ageing_fixed!(population, pi, N)

    new_pi = zeros(population.m,population.ageStructure.Amax)
    new_N = zeros(population.ageStructure.Amax)
    for a in 2:population.ageStructure.Amax
        
        
        pi_a = population.pi[:,a-1]
        Npi_a = population.abundance[a-1] .* pi_a
        sa = population.ageStructure.Survival[a-1]
       
        Npi_a = broadcast(n -> sampleN(round(Int,n),sa),Npi_a)
        
        if sum(Npi_a) > 0
            new_pi[:,a] = Npi_a / sum(Npi_a)
        else
            new_pi[:,a] = zeros(population.m)
                
        end 
        new_N[a] = sum(Npi_a)
    end 
    
    new_pi[:,1] = pi
    new_N[1] = N 
 
    population.abundance = new_N
    population.pi= new_pi
    
end 




function time_step_fixed_DSI!(population, N_im)
    pi, N = reproduction_fixed(population)
    pi, N =  recruitment_fixed(pi, N, population, true)
    pi, N = selection_fixed(pi,N,population)
    pi, N = immigration_fixed(pi,N,N_im)
    ageing_fixed!(population,pi, N)
end 



function fitness_fixed(population)
    
    pi,N = reproduction_fixed(population)

    z = population.z
    theta = population.theta
    s = population.s
    V = population.V
    V_prime = (s+1/V)^(-1)
    z_prime = V_prime.*(s*theta .+ z./V)
    p_survival = broadcast(i -> prop_survival(z[i],z_prime[i],s,theta,V, V_prime),1:length(z))

    return sum(p_survival .* pi)
    
end 


end 