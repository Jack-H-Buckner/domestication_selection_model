"""
calculates fecundity at age relationships based on
Von-bertalanffy growth curve
"""
module fecundity_at_age

function isometric_LVB(a,f_infty,k,a0)
    F = f_infty*(1-exp(-k*(a - a0)))^3
    if F > 0
        return F
    else
        return 0
    end
end 

allometric_LVB(a,f_infty,k,a0,beta) = f_infty*(1-exp(-k*(a - a0)))^beta



function isometric_LVB_A_mature(a,A_mature,f_infty,k,a0)
    v = f_infty*(1-exp(-k*(a - a0)))^3
    if a > A_mature &&  v>0 
        return f_infty*(1-exp(-k*(a - a0)))^3
    else
        return 0
    end
end 

end #module 