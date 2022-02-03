"""
Base parameter set for coupled genetic demographic model 

uses parameters from Smyth 2016 for demographic processes and quantitative trait to represent fittness effects

"""
module coupled_genetic_model
include("../age_structure/full_trait_distribution.jl")
include("../White_sturgeon/demo_params.jl")
using Roots
using Distributions

## model1 
s1 = 0.1
theta = 0.0
Vle = 1.0
k_ind = 3
R_star  = demo_params.Smyth_16_R_star
alpha = demo_params.Smyth_2016_alpha[k_ind]
beta = demo_params.Smyth_2016_beta[k_ind]
A_max = demo_params.A_max
survival = demo_params.Smyth_2016_survival
fecundity = demo_params.Smyth_16_F

model1 = full_trait_distribution.init_population_2(A_max,survival,fecundity, R_star, 
                                                    alpha , beta, theta, 
                                                    s1, -5, 15, 0.05, Vle, true)


theta_1 = 5.0
d = Distributions.Normal(theta_1,sqrt(Vle))
dsn_im = pdf.(d,model1.grid)
dsn_im_1 = dsn_im ./ sum(dsn_im)

function calc_RRS(pop, theta, Vle)
    W0 = sum(pop.trait[:,1].*pop.gradient)
    dsn1 = Distributions.Normal(theta,sqrt(Vle))
    dsn1 = pdf.(dsn1,pop.grid)
    dsn1 = dsn1 ./ sum(dsn1)
    W1 = sum(dsn1.*pop.gradient)
    return W1/W0
end 

function init_model1(s,RRS, k_ind, correction)
    # init populaiton 
    alpha = demo_params.Smyth_2016_alpha[k_ind]
    beta = demo_params.Smyth_2016_beta[k_ind]
    pop = full_trait_distribution.init_population_2(A_max,survival,fecundity, R_star, 
                                                    alpha , beta, theta, 
                                                    s, -5, 15, 0.05, Vle,correction)
    # solve for migrants DSN
    f = x -> calc_RRS(pop, x, Vle) - RRS
    theta_1 = Roots.find_zero(f, (0.0, 12.0))
    d = Distributions.Normal(theta_1,sqrt(Vle))
    dsn_im = pdf.(d,pop.grid)
    dsn_im = dsn_im ./ sum(dsn_im)
    
    return pop, dsn_im, floor(Int,demo_params.Smyth_2016_T1 )
end 





"""
selection -> imigraiton -> density dependence
"""
function update_SID!(pop, R_im, dsn_im)

    # reproduction
    dsn, F = full_trait_distribution.reproduction(pop)
    # seleciton
    dsn, R = full_trait_distribution.selection(dsn, F, pop)
    # immigration
    dsn, R = full_trait_distribution.immigration(dsn, R, dsn_im, R_im)
    # dentisy dependence
    R = full_trait_distribution.recruitment(R, pop)
    # ageing!
    pop = full_trait_distribution.ageing!(pop, R, dsn)
    return pop, R
end


"""
selection -> density dependence -> imigraiton
"""
function update_SDI!(pop, R_im, dsn_im)

    # reproduction
    dsn, F = full_trait_distribution.reproduction(pop)
    # seleciton
    dsn, R = full_trait_distribution.selection(dsn, F, pop)
    # dentisy dependence
    R = full_trait_distribution.recruitment(R, pop)
    # immigration
    dsn, R = full_trait_distribution.immigration(dsn, R, dsn_im, R_im)
    # ageing!
    pop = full_trait_distribution.ageing!(pop, R, dsn)
    return pop, R
end

"""
imigraiton -> selection -> density dependence
"""
function update_ISD!(pop, R_im, dsn_im)

    # reproduction
    dsn, F = full_trait_distribution.reproduction(pop)
    # immigration
    dsn, R = full_trait_distribution.immigration(dsn, F, dsn_im, R_im)
    # seleciton
    dsn, R = full_trait_distribution.selection(dsn, R, pop)
    # dentisy dependence
    R = full_trait_distribution.recruitment(R, pop)
    # ageing!
    pop = full_trait_distribution.ageing!(pop, R, dsn)
    return pop, R
end

"""
imigraiton -> density dependence -> selection 
"""
function update_IDS!(pop, R_im, dsn_im)

    # reproduction
    dsn, R = full_trait_distribution.reproduction(pop)
    # immigration
    dsn, R = full_trait_distribution.immigration(dsn, R, dsn_im, R_im)
    # dentisy dependence
    R = full_trait_distribution.recruitment(R, pop)
    # seleciton
    dsn, R = full_trait_distribution.selection(dsn, R, pop)
    # ageing!
    pop = full_trait_distribution.ageing!(pop, R, dsn)
    return pop, R
end


"""
imigraiton -> density dependence -> selection 
"""
function update_DSI!(pop, R_im, dsn_im)

    # reproduction
    dsn, R = full_trait_distribution.reproduction(pop)
    # dentisy dependence
    R = full_trait_distribution.recruitment(R, pop)
    # seleciton
    dsn, R = full_trait_distribution.selection(dsn, R, pop)

    # immigration
    dsn, R = full_trait_distribution.immigration(dsn, R, dsn_im, R_im)
    # ageing!
    pop = full_trait_distribution.ageing!(pop, R, dsn)
    return pop, R
end
end # module 