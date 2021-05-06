import eos
methane = eos.which_mol("CH4",190.4 , 46.0, 0.011)
ethane = eos.which_mol("C2H6", 305.4, 48.8, 0.099)

T = 300 # K 
P_total = 1.01325 # bar 
x = [0.2, 0.8] # mole fractions
sigma = - 0.004 # binary interaction parameter for Xe/Kr

sol = eos.eosPengR_mix(methane, ethane, T, P_total, x, sigma)
