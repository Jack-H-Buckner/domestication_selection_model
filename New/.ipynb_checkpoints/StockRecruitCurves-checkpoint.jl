module StockRecruitCurves

using Roots

mutable struct BevetonHolt
    a::Float64
    b::Float64
    R_star::Float64
    k::Float64
    LEP::Float64
end


function (sr::BevetonHolt)(x)
    return sr.a*x/(1+sr.b*x)
end 

"""
    init_BevetonHolt(R_hat,k)

initalizes beverton-holt curve with compensaton ration k and
equilibrium recruitment R_hat, givin 
"""
function init_BevetonHolt(R_star::Float64,k::Float64, LEP::Float64)
    a = k ./ LEP
    b = (k.-1) ./ (R_star * LEP)
    BevetonHolt(a,b,R_star,k,LEP)
end 

function update_BevetonHolt_k!(sr, k::Float64)
    sr.a = k ./ sr.LEP
    sr.b = (k.-1) ./ (sr.R_star * sr.LEP)
    sr.k = k
end 

function update_BevetonHolt_Rstar!(sr, R_star::Float64)
    sr.a = sr.k ./ sr.LEP
    sr.b = (sr.k.-1) ./ (R_star * sr.LEP)
    sr.R_star = R_star
end 

function update_BevetonHolt_LEP!(sr, LEP::Float64)
    sr.a = sr.k ./ LEP
    sr.b = (sr.k.-1) ./ (R_star * LEP)
    sr.LEP = LEP
end 

"""
    quilibrium_size(survival::AbstractVector{Float64}, R_star::Float64, Phi::Function )

survival - vertor probility of surviving from age a to a+1
R_star - Float recruitment rate
Phi - Function if supplied weights the population by age
"""
function equilibrium_size(survival::AbstractVector{Float64}, R_star::Float64, Phi::Function )
    N = R_star
    acc = 0
    for i in 1:length(survival)
        acc += Phi(i)*N
        N *= survival[i]*N
    end 
    return acc
end 

function equilibrium_size(survival::AbstractVector{Float64}, R_star::Float64)
    N = R_star
    acc = 0
    for i in 1:length(survival)
        acc += N
        N *= survival[i]*N
    end 
    return acc
end 

"""
    solve_Rstar(N::Float64, survival::AbstractVector{Float64})

N - float equilibrium population size 
N - survival rate
Amat - int if provided is used as age of maturation, and will assume N is tartegt size of matru population
Phi - funtion if provided is used as a age specific weight for calcualting pop size 
"""
function solve_Rstar(N::Float64, survival::AbstractVector{Float64})
    f = x -> equilibrium_size(survival, x) - N
    R_star = Roots.find_zero(f, 0, 100.0*N)
    return R_star
end 

function solve_Rstar(N::Float64, survival::AbstractVector{Float64}, Amat::Float64)
    Phi = x -> 1. * (x<Amat)
    f = x -> equilibrium_size(survival, x, Phi) - N
    R_star = Roots.find_zero(f, 0, 100.0*N)
    return R_star
end 

function solve_Rstar(N::Float64, survival::AbstractVector{Float64}, Phi::Function)
    f = x -> equilibrium_size(survival, x, Phi) - N
    R_star = Roots.find_zero(f, 0, 100.0*N)
    return R_star
end 



end 