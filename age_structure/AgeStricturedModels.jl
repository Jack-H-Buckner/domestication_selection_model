"""
This module contains the basic machinery required to analyze
the coupled age structured genetic demogrpahic model. For now,
the model tracks the number of individuals at a given age and the 
mean trait value for the cohort. 
"""
module AgeStricturedModels


"""
This struct stores all of the paramters needed to 
define transition matrices for the abudnaces and 
mean trait value of a cohort. 
The model assume that selection effects survival between 
age classes.  
"""
mutable struct AgeStricturedModel    
    A_max::Int
    N::AbstractVector{Float64}
    z::AbstractVector{Float64}
    fecundity::AbstractVector{Float64}
    survival::AbstractVector{Float64}
    selection::AbstractVector{Float64}
end 

function build_age_matrix(A_max,fecundity,survival)
    A = zeros(A_max,A_max)
    A[1,:] = fecundity
    for i in 1:(A_max-1)
        A[i+1,i] = survival[i]
    end 
    return A
end 

function build_matrices_linear(LMS)
    N_mat = zeros(LMS.A_max,LMS.A_max)
    z_mat = zeros(LMS.A_max,LMS.A_max)
    ## survival 
    survival = LMS.survival.*exp.(-1 .*LMS.z.*LMS.selection./2)
    selection = 1 ./(1 .+ LMS.selection)
    ## new triats 
    F = sum(LMS.N.*LMS.fecundity)
    z_weights = LMS.N.*LMS.fecundity./F
    ## fill in matrices
    N_mat[1,:] = LMS.fecundity
    z_mat[1,:] = z_weights
    for i in 1:(LMS.A_max-1)
        N_mat[i+1,i] = survival[i]
        z_mat[i+1,i] = selection[i]
    end 
    return N_mat, z_mat
end


function time_step_linear!(LMS)
    N_mat,z_mat = build_matrices_linear(LMS)
    LMS.z = z_mat*LMS.z
    LMS.N = N_mat*LMS.N
    
    return LMS
end


function total_fecundity(LMS)
    F = sum(LMS.N.*LMS.fecundity)
    return F
end 


function percapita_fecundity(LMS)
    F = sum(LMS.N.*LMS.fecundity)./sum(LMS.N)
    return F
end 

function build_matrices_afine(LMS)
    N_mat = zeros(LMS.A_max,LMS.A_max)
    z_mat = zeros(LMS.A_max,LMS.A_max)
    ## survival 
    survival = LMS.survival.*exp.(-1 .*LMS.z.*LMS.selection./2)
    selection = 1 ./(1 .+ LMS.selection)
    ## new triats 
    F = sum(LMS.N.*LMS.fecundity)
    z_weights = LMS.N.*LMS.fecundity./F
    ## fill in matrices
    z_mat[1,:] = z_weights
    for i in 1:(LMS.A_max-1)
        N_mat[i+1,i] = survival[i]
        z_mat[i+1,i] = selection[i]
    end 
    return N_mat, z_mat
end


function time_step_afine!(LMS, R)
    N_mat,z_mat = build_matrices_afine(LMS)
    LMS.z = z_mat*LMS.z
    LMS.N = N_mat*LMS.N .+ R
    
    return LMS
end


function age_distribution_given_recruitment(R, survival, A_max )
    N = zeros(A_max)
    N[1] = R
    for i in 2:A_max
        N[i] = N[i-1]*survival[i-1]
    end 
    return N
end 

function LEP(R, survival, fecundity, A_max )
    N = age_distribution_given_recruitment(R, survival, A_max )
    return sum(N.*fecundity)
end 



end # module 