"""
calculates fecundity at age relationships based on
Von-bertalanffy growth curve
"""
module fecundity_at_age

isometric_LVB(a,f_infty,k,a0) = f_infty*(1-exp(-k*(a - a0)))^3
allometric_LVB(a,f_infty,k,a0,beta) = f_infty*(1-exp(-k*(a - a0)))^beta



function isometric_LVB_A_mature(a,A_mature,f_infty,k,a0)
    if a > A_mature
        return f_infty*(1-exp(-k*(a - a0)))^3
    else
        return 0
    end
end 

end #module 