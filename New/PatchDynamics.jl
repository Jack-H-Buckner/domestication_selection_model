"""
    PatchDynamics

Combines populaitons from `AgeTraitStructuredModels.jl` into a network of patches to represent discrete
spatial structure with gene flow between patchs
"""
module PatchDynamics

include("AgeTraitStructuredModels.jl")

mutable struct Network{T1}
    populations
    juvinile_dispersal
    adult_dispersal
    reproduction
end 




end # module 