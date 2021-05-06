import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

#computes the Peng-Robinson equation of state
R = 8.314e-5 #universal gas const in m3-bar/K-mol

class which_mol:

	#specifying the molecule
	def __init__(self, name, Tc, Pc, w):
		self.name = name #molecule name
		self.Tc = Tc #critical temperature (K)
		self.Pc = Pc #critical pressure (bar)
		self.w = w #accentric factor

	#printing parameters
	def print_params(self):
		print("Mol. %s, Tc = %.1f K, Pc = %.1f bar, w =%f" % (self.name, self.Tc, self.Pc, self.w))

def eosPengR(which_mol, T, P, plot=True, results=True):
	#parameterts in terms of critical constants and acentric factor
	Tr = T / which_mol.Tc #reduced temperature
	a = 0.45724 * R * R * which_mol.Tc * which_mol.Tc / which_mol.Pc
	b = 0.07780 * R * which_mol.Tc / which_mol.Pc
	k = 0.37464 + 1.54226 * which_mol.w - 0.26992 * which_mol.w * which_mol.w
	alpha = (1 + k * (1 - np.sqrt(Tr)))**2.0
	A = alpha * a * P / R**2.0 / T**2.0
	B = b * P / R / T

	#the polynomial we need to solve to get the compressibility
	def fn(Z):
		return Z**3.0 - (1-B) * Z**2.0 + (A-2*B-3*B**2) * Z - (A*B - B*B - B**3)
	Z = newton(fn, 0.9) #compressibility

	phi = np.exp(Z-1-np.log(Z-B)-A/np.sqrt(8.0)/B * np.log((Z+(np.sqrt(2.0)+1)*B)/(Z-(np.sqrt(2.0)-1)*B))) #fugacity ceofficient
	rho = P /R/T/Z #density

	#print results if "True"
	if results:
		print("Caculations done at T = %.1f K, P = %.1f bar" %(T, P))
		print("Compressibility Z = %.9f, desnity, rho = %.9f mol/m3, fugacity coefficient phi = %.9f" %(Z, rho, phi))

	#plot if "True"
	if plot:
		x = np.linspace(-0.5, 1.5, 50)
		fig = plt.figure()
		plt.plot(x, fn(x), color = 'r')
		plt.xlabel('Compressibility Z')
		plt.ylabel('fn(Z)')
		plt.axvline(x=Z)
		plt.axhline(y=0)
		plt.show()	

#Peng-Robinson EOS for a binary mixture
def eosPengR_mix(which_mol_a, which_mol_b, T, P, x, sigma, plot = True, results = True):

	#array containing Tc and Pc of the two components

	Tc = np.array([which_mol_a.Tc, which_mol_b.Tc])
	Pc = np.array([which_mol_a.Pc, which_mol_b.Pc])
	W = np.array([which_mol_a.w, which_mol_b.w])

	#x denotes the mole fraction
	#am and bm are the mixing parameters
	#sigma is the interaction parameter

	Tr = T / Tc
	K = 0.37464 + 1.54226 * W - 0.26992 * W * W
	tempa = 0.457235 * R**2 * Tc**2 / Pc
	a = tempa * (1 + K * (1-np.sqrt(Tr)))**2.0
	b = 0.07780 * R * Tc / Pc
	aij = (1-sigma) * np.sqrt(a[0] * a[1])
	am = a[0] * x[0]**2 + 2 * aij * x[0] * x[1] + a[1] * x[1]**2
	bm = x[0] * b[0] + x[1] * b[1]
	A = am * P/R/R/T/T
	B = bm * P/R/T

	def fn(Z):
		return Z**3.0 - (1-B) * Z**2.0 + (A-2*B-3*B**2) * Z - (A*B - B*B - B**3)
	Z = newton(fn, 0.9)
	
	lphi_0 = -np.log(Z - B) + (Z - 1.0) * b[0] / bm - A / np.sqrt(8) / B * (2.0 / am * (x[0] * a[0] + x[1] * aij) - b[0] / bm) * np.log((Z + (1.0 + np.sqrt(2)) * B) / (Z + (1.0 - np.sqrt(2)) * B))
	lphi_1 = -np.log(Z - B) + (Z - 1.0) * b[1] / bm - A / np.sqrt(8) / B * (2.0 / am * (x[1] * a[1] + x[0] * aij) - b[1] / bm) * np.log((Z + (1.0 + np.sqrt(2)) * B) / (Z + (1.0 - np.sqrt(2)) * B))

	PHI = np.exp(np.array([lphi_0, lphi_1]))

	#print results if "True"
        if results:
                print("Caculations done at T = %.1f K, P = %.1f bar" %(T, P))
		print("Com01: %s fugacity coefficient phi = %f" %(which_mol_a.name, PHI[0]))
                print("Com02: %s fugacity coefficient phi = %f" %(which_mol_b.name, PHI[1]))

        #plot if "True"
        if plot:
                x = np.linspace(-0.5, 1.5, 50)
                fig = plt.figure()
                plt.plot(x, fn(x), color = 'r')
                plt.xlabel('Compressibility Z')
                plt.ylabel('fn(Z)')
                plt.axvline(x=Z)
                plt.axhline(y=0)
                plt.show()  
