# opticalmodel_deuterons
# Functions for deuteron optical models
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# =============================================================================================== #
from opticalmodel_globals import *
import numpy as np

# Define optical model dictionary
deuteron_dct = {
	"AC":  [0,1],
	"B":   [1,2],
	"DNR": [2,3],
	"DR":  [3,4],
	"HSS": [4,5],
	"LH":  [5,6],
	"PP":  [6,7],
	"len":     7,
	"ALL-D":[0,7]
}

# Return the names of the deuteron potentials alphabetically
def DeuteronModelNumber():
	return [ "AnCai", "Bojowald", "DaehnickNR", "DaehnickR", "HanShiShen", "LohrHaeberli", "PereyPerey" ] 

# An and Cai
def AnCai(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	v = 91.85 - 0.249*E + (1.16e-4)*(E**2) + 0.642*Z*(A**(-1.0/3.0))
	vi = 1.104 + 0.0622*E
	vsi = 10.83 - 0.0306*E
	vso = 3.557
	vsoi = 0
	
	r0 = 1.152 - 0.00776*(A**(-1.0/3.0))
	ri0 = 1.305 + 0.0997*(A**(-1.0/3.0))
	rsi0 = 1.334 + 0.152*(A**(-1.0/3.0))
	rso0 = 0.972
	rsoi0 = 0
	
	a = 0.719 + 0.0126*(A**(1.0/3.0))
	ai = 0.855 - 0.1*(A**(1.0/3.0))
	asi = 0.531 + 0.062*(A**(1.0/3.0))
	aso = 1.011
	asoi = 0
	
	rc0 = 1.303

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "An and Cai deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)
	return string_list

# =============================================================================================== #
# PereyPerey is another potential used for deuterons
def PereyPerey(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
	
	# CALCULATE FINAL PARAMETERS
	v = 81 - 0.22*E + 2*Z*(A**(-1.0/3.0))
	vi = 0.0
	vsi = 14.4 + 0.24*E
	vso = 0.0
	vsoi = 0.0
	
	r0 = 1.15
	ri0 = 0.0
	rsi0 = 1.34
	rso0 = 0.0
	rsoi0 = 0.0
	
	a = 0.81
	ai = 0.0
	asi = 0.68
	aso = 0.0
	asoi = 0.0
	
	rc0 = 1.15

	# Format final paramaters into a list of strings
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Perey and Perey deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)
	return string_list

# =============================================================================================== #
# LohrHaeberli is another potential used for deuterons
def LohrHaeberli(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	v = 91.13 + 2.2*Z*(A**(-1.0/3.0))
	vi = 0.0
	vsi = 218*(A**(-2.0/3.0))
	vso = 7.0
	vsoi = 0.0
	
	r0 = 1.05
	ri0 = 0.0
	rsi0 = 1.43
	rso0 = 0.75
	rsoi0 = 0.0
	
	a = 0.86
	ai = 0.0
	asi = 0.5 + 0.013*(A**(2.0/3.0))
	aso = 0.5
	asoi = 0.0
	
	rc0 = 1.30

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Lohr and Haeberli deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)

	return string_list

# =============================================================================================== #
# HanShiShen is another potential used for deuterons
def HanShiShen(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	v = 82.18 - 0.148*E - 0.000886*E*E - 34.811*( N - Z )/A + 1.058*Z*( A**(-1.0/3.0) )
	vi = -4.916 + 0.0555*E + 0.0000442*E*E + 35.0*( N - Z )/A if -4.916 + 0.0555*E + 0.0000442*E*E + 35.0*( N - Z )/A > 0.0 else 0.0
	vsi = 20.968 - 0.0794*E - 43.398*( N - Z )/A
	vso = 3.703
	vsoi = -0.206
	
	r0 = 1.174
	ri0 = 1.563
	rsi0 = 1.328
	rso0 = 1.234
	rsoi0 = 1.234
	
	a = 0.809
	ai = 0.7 + 0.045*( A**(1.0/3.0) )
	asi = 0.465 + 0.045*( A**(1.0/3.0) )
	aso = 0.813
	asoi = 0.813
	
	rc0 = 1.698

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Han, Shi, and Shen deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)

	return string_list

# =============================================================================================== #
# Bojowald is another potential used for deuterons
def Bojowald(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	v = 81.33 + 1.43*Z*( A**(-1.0/3.0) ) - 0.24*E
	vi = 0.132*( E - 45.0 ) if 0.132*( E - 45.0 ) > 0.0 else 0.0
	vsi = 7.8 + 1.04*( A**(1.0/3.0) ) - 0.712*vi
	vso = 6.0
	vsoi = 0.0
	
	r0 = 1.18
	ri0 = 1.27
	rsi0 = 1.27
	rso0 = 0.78 + 0.038*( A**(1.0/3.0) )
	rsoi0 = 0.0
	
	a = 0.636 + 0.035*( A**(1.0/3.0) )
	ai = 0.768 + 0.021*( A**(1.0/3.0) )
	asi = 0.768 + 0.021*( A**(1.0/3.0) )
	aso = 0.78 + 0.038*( A**(1.0/3.0) )
	asoi = 0.0
	
	rc0 = 1.3

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Bojowald deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)

	return string_list

# =============================================================================================== #
# DaehnickNR is another potential used for deuterons
def DaehnickNR(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	mu = [ 8, 20, 28, 50, 82, 126 ]
	step1 = [ ( 0.5*(ii - N) )**2 for ii in mu ]
	step2 = [ np.exp(-jj) for jj in step1 ]
	step3 = sum(step2)
	beta = -1.0*( ( 0.01*E )**2 )

	v = 88.5 - 0.26*E + 0.88*Z*( A**(-1.0/3.0) )
	vi = ( 12.2 + 0.026*E )*( 1.0 - np.exp(beta) )
	vsi = ( 12.2 + 0.026*E )*np.exp(beta)
	vso = 7.33 - 0.029*E
	vsoi = 0.0
	
	r0 = 1.17
	ri0 = 1.325
	rsi0 = 1.325
	rso0 = 1.07
	rsoi0 = 0.0
	
	a = 0.709 + 0.0017*E
	ai = 0.53 + 0.07*( A**(1.0/3.0) ) - 0.04*step3
	asi = 0.53 + 0.07*( A**(1.0/3.0) ) - 0.04*step3
	aso = 0.66
	asoi = 0.0
	
	rc0 = 1.3

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Daehnick (non-rel) deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)
		print( "beta\t" + str(beta) )
		print( "mu\t" + str(mu) )
		print( "step1\t" + str(step1) )
		print( "step2\t" + str(step2) )
		print( "step3\t" + str(step3) )
	return string_list

# =============================================================================================== #
# DaehnickR is another potential used for deuterons
def DaehnickR(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
		
	# CALCULATE FINAL PARAMETERS
	mu = [ 8, 20, 28, 50, 82, 126 ]
	step1 = [ ( 0.5*(ii - N) )**2 for ii in mu ]
	step2 = [ np.exp(-jj) for jj in step1 ]
	step3 = sum(step2)
	beta = -1.0*( ( 0.01*E )**2 )

	v =  88.0 - 0.283*E + 0.88*Z*( A**(-1.0/3.0) )
	vi = ( 12 + 0.031*E )*( 1 - np.exp(beta) )
	vsi = ( 12 + 0.031*E )*np.exp(beta)
	vso = 7.2 - 0.032*E
	vsoi = 0.0
	
	r0 = 1.17
	ri0 = 1.376 - 0.01*np.sqrt(E)
	rsi0 = 1.376 - 0.01*np.sqrt(E)
	rso0 = 1.07
	rsoi0 = 0.0
	
	a = 0.717 + 0.0012*E
	ai = 0.52 + 0.07*( A**(1.0/3.0) ) - 0.04*step3
	asi = 0.52 + 0.07*( A**(1.0/3.0) ) - 0.04*step3
	aso = 0.66
	asoi = 0.0
	
	rc0 = 1.3

	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Daehnick (non-rel) deuteron")
		PrintCalculatedQuantities(A,Z,E,Q)
		print( "beta\t" + str(beta) )
		print( "mu\t" + str(mu) )
		print( "step1\t" + str(step1) )
		print( "step2\t" + str(step2) )
		print( "step3\t" + str(step3) )
	return string_list

# 124Te(p,d), Ebeam=22MeV, Ex=0MeV
#PRINT = 1
#DaehnickR(123, 52, 22, 0, 123.9028179, 1.00782503224, 2.01410177811, 122.9042698, 1)
