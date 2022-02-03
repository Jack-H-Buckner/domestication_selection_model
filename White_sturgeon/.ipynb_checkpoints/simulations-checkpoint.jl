"""
code for simualtions 
"""
module simulations

include("coupled_genetic_model.jl")
include("../age_structure/full_trait_distribution.jl")
function simulation_DSI(R_im, duratation, N_after, s, RRS, k_ind)
    
    # initialize populations 
    pop, dsn_im, T_ = coupled_genetic_model.init_model1(s,RRS, k_ind, true)
    
    N = floor(Int,T_ + duratation*T_ + N_after * T_)
    W = zeros(N)
    E_ = zeros(N)
    R_im = vcat(repeat([0], T_), repeat([R_im], floor(Int,duratation*T_)), repeat([0.0], floor(Int,N_after*T_)))

    for t in 1:N
        pop, R = coupled_genetic_model.update_DSI!(pop, R_im[t], dsn_im)
        W[t] = full_trait_distribution.fittness(pop)
        E_[t] = sum(pop.abundance.*pop.fecundity)
    end 
    return W, E_, R_im, T_
end


function quantities(W,E,R_im, T_)
    min_fittness = W[argmax(-1 .* W)]
    min_abundance = E[argmax(-1 .* E)]
    T_min_abundnace = argmax(-1 .* E) * T_
    integrated_fitness = sum(W .-1) ./ (length(R_im) .- T_)
    return min_fittness, min_abundance, T_min_abundnace, integrated_fitness
end 


end # module 