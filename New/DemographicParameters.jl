module DemographicParameters

include("FecundityAgeRelationships.jl")
include("StockRecruitCurves.jl")
##########################
#### fecundity at age ####
##########################

A_max = 150
a = collect(1:A_max)

# Wood et al 2007
Wood_07_W100 = 7
Wood_07_Wmat = 165^3 *Wood_07_W100/100
Wood_07_K = 0.023
Wood_07_Linfty = 412
Wood_07_mu_s1 = 80# 100 # 
Wood_07_sigma_s1 = 5 #0 #5
Wood_07_mu_s2 = 200  
Wood_07_sigma_s2 = 0 #0 #5
Wood_07_m = 0.923


# Smyth 2016
Smyth_16_W100 = 9.3
Smyth_16_Wmat = 153^3 *Smyth_16_W100/100
Smyth_16_K = 0.034
Smyth_16_Linfty = 300
Smyth_16_mu_s = 200  
Smyth_16_sigma_s = 0 #0 #5
Smyth_16_m = 0.96


# fecundity age function
Wood_07_F1 = FecundityAgeRelationships.Walters2006Senescence(Wood_07_Linfty, Wood_07_W100, Wood_07_Wmat,
                                                            Wood_07_K, 0, Wood_07_mu_s1, Wood_07_sigma_s1)
Wood_07_F1_vec = Wood_07_F1.(a)

Wood_07_F2 = FecundityAgeRelationships.Walters2006(Wood_07_Linfty, Wood_07_W100, Wood_07_Wmat,Wood_07_K, 0)
Wood_07_F2_vec = Wood_07_F2.(a)

Smyth_16_F = FecundityAgeRelationships.Walters2006(Smyth_16_Linfty, Smyth_16_W100, Smyth_16_Wmat, Smyth_16_K, 0)
Smyth_16_F_vec = Smyth_16_F.(a)

##########################
####  survival rates  ####
##########################

Wood_2007_survival = repeat([Wood_07_m],A_max)
Smyth_2016_survival = repeat([Smyth_16_m],A_max) 
long_lived_survival = vcat(repeat([0.9],80), repeat([0.98],A_max - 80))
    
############################
####  impact functions  ####
############################


# life time egg production 
Wood_2007_LEP1 = sum(Wood_07_F1_vec .* Wood_07_m .^a)
Wood_2007_impact1 = Wood_07_F1_vec  .* Wood_07_m .^a

Wood_2007_LEP2 = sum(Wood_07_F2_vec  .* Wood_07_m .^a)
Wood_2007_impact2 = Wood_07_F2_vec  .* Wood_07_m .^a

Smyth_2016_LEP = sum(Smyth_16_F_vec  .* Smyth_16_m .^a)
Smyth_2016_impact = Smyth_16_F_vec  .* Smyth_16_m .^a

psurvival = zeros(A_max)
psurvival[1] = long_lived_survival[1]
for i in 2:A_max
    psurvival[i] = long_lived_survival[i] * psurvival[i-1]
end 

Wood_2007_LEP1_ll = sum(Wood_07_F1_vec  .* psurvival)
Wood_2007_impact1_ll = Wood_07_F1_vec  .* psurvival

Wood_2007_LEP2_ll = sum(Wood_07_F2_vec  .* psurvival)
Wood_2007_impact2_ll = Wood_07_F2_vec  .* psurvival

Smyth_2016_LEP_ll = sum(Smyth_16_F_vec  .* psurvival)
Smyth_2016_impact_ll = Smyth_16_F_vec  .* psurvival





##########################
####  stock recruit   ####
##########################



N_star = 1500 # equilibirum popuatlion size 
Smyth_16_R_star = N_star/sum(Smyth_16_m .^a) # equilibrium recruitment rate
Wood_2007_R_star = N_star/sum(Wood_07_m .^a)
LL_R_star = N_star/sum(psurvival)
k = 5.0 #[1.5,3,5,7,9]

Smyth_2016_sr = StockRecruitCurves.init_BevetonHolt(Smyth_16_R_star,k, Smyth_2016_LEP)
Smyth_2016_ll_sr = StockRecruitCurves.init_BevetonHolt(LL_R_star,k, Smyth_2016_LEP_ll)
Wood_2007_1_sr = StockRecruitCurves.init_BevetonHolt(Wood_2007_R_star,k, Wood_2007_LEP1)
Wood_2007_2_sr = StockRecruitCurves.init_BevetonHolt(Wood_2007_R_star,k, Wood_2007_LEP2)
Wood_2007_1_ll_sr = StockRecruitCurves.init_BevetonHolt(LL_R_star,k, Wood_2007_LEP1_ll)
Wood_2007_2_ll_sr = StockRecruitCurves.init_BevetonHolt(LL_R_star,k, Wood_2007_LEP2_ll)



############################
####  generation time   ####
############################

Wood_2007_T1 =  sum(a .* (Wood_07_F1_vec  .* Wood_07_m .^a) ./ Wood_2007_LEP1)
Wood_2007_T2 =  sum(a .* (Wood_07_F2_vec  .* Wood_07_m .^a) ./ Wood_2007_LEP2)
Smyth_2016_T1 =  sum(a .* (Smyth_16_F_vec  .* Smyth_16_m .^a) ./ Smyth_2016_LEP)

Wood_2007_T3 =  sum(a .* (Wood_07_F1_vec   .* psurvival) ./ Wood_2007_LEP1_ll)
Wood_2007_T4 =  sum(a .* (Wood_07_F2_vec   .* psurvival) ./ Wood_2007_LEP2_ll)
Smyth_2016_T2 =  sum(a .* (Smyth_16_F_vec  .* psurvival) ./ Smyth_2016_LEP_ll)

end # modul 