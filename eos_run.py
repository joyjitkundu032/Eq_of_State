import eos
#Critical Constants and Acentric Factors for Selected Fluids
C6H6 = eos.which_mol("Benzene", 562.1, 48.9, 0.212)
C6H6.print_params()

#Solving for the compressibility factor 
sol = eos.eosPengR(C6H6, 300.0, 10.0, plot=True, results=True)
