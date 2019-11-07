# opticalmodel_protons
# Functions for proton optical models
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# =============================================================================================== #
from opticalmodel_globals import *
import numpy as np

# Define optical model dictionary
proton_dct = {
	"BG":    [0,1],
	"KD":    [1,2],
	"M":     [2,3],
	"P":     [3,4],
	"V":     [4,5],
	"len":       5,
	"ALL-P": [0,5]
}

# Return the names of the proton potentials alphabetically
def ProtonModelNumber():
	return [ "BecchettiGreenlees", "KoningDelaroche", "Menet", "Perey", "Varner" ]

# =============================================================================================== #
# Koning-Delaroche set
def KoningDelaroche(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
	
		
	# Calculate more involved parameters
	vp1 = 59.3 + 21*float((N-Z))/float(A) - 0.024*A
	vp2 = 0.007067 + (4.23e-6)*A
	vp3 =  (1.729e-5) + (1.136e-8)*A
	vp4 = 7e-9
	
	wp1 = 14.667 + 0.009629*A
	wp2 = 73.55 + 0.0795*A
	
	dp1 = 16*(1 + float(N-Z)/float(A))
	dp2 = 0.018 + 0.003802/( 1 + np.exp( (A - 156)/8 ) )
	dp3 = 11.5
	
	vpso1 = 5.922 + 0.0030*A
	vpso2 = 0.0040
	
	wpso1 = -3.1
	wpso2 = 160
	
	epf = -8.4075 + 0.01378*A
	rc = 1.198 + 0.697*(A**(-2.0/3.0)) + 12.994*(A**(-5.0/3.0))
	
	vc = 1.73*Z*(A**(-1.0/3.0))/rc
	
	# Calculate final parameters
	v = vp1*( 1 - (vp2*(E - epf)) + (vp3*((E-epf)**2)) - (vp4*((E-epf)**3)) ) + ( vc*vp1*( vp2 - (2*vp3*(E-epf)) + (3*vp4*((E-epf)**2)) ) )
	vi = wp1*((E-epf)**2)/(((E-epf)**2) + (wp2**2))
	vsi = dp1*((E-epf)**2)/(((E-epf)**2) + (dp3**2))*np.exp( -dp2*(E-epf) )
	vso = vpso1*np.exp( -vpso2*(E-epf) )
	vsoi = wpso1*((E-epf)**2)/(((E-epf)**2) + (wpso2**2))
	
	r0 = 1.3039 - 0.4054*(A**(-1.0/3.0))
	ri0 = r0
	rsi0 = 1.3424 - 0.01585*(A**(1.0/3.0))
	rso0 = 1.1854 - 0.647*(A**(-1.0/3.0))
	rsoi0 = rso0
	
	a = 0.6778 - 0.0001487*A
	ai = a
	asi = 0.5187 + 0.0005205*A
	aso = 0.59
	asoi = aso
	
	rc0 = rc
	
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Koning and Delaroche ")
		PrintCalculatedQuantities(A,Z,E,Q)
		print(DIV)
		print("vp1" + sep + str(vp1))
		print("vp2" + sep+ str(vp2))
		print("vp3" + sep + str(vp3))
		print("vp4" + sep + str(vp4))
		print(div)
		print("wp1" + sep + str(wp1))
		print("wp2" + sep + str(wp2))
		print("dp1" + sep + str(dp1))
		print("dp2" + sep + str(dp2))
		print("dp3" + sep + str(dp3))
		print(div)
		print("vpso1" + sep + str(vpso1))
		print("vpso2" + sep + str(vpso2))
		print("wpso1" + sep + str(wpso1))
		print("wpso2" + sep + str(wpso2))
		print(div)
		print("epf" + sep + str(epf))
		print("rc" + sep + str(rc))
		print("vc" + sep + str(vc))
		print(DIV)

	return string_list

# =============================================================================================== #
# Perey (protons)
def Perey(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
	
	# Calculate final parameters
	v = 53.3 - 0.55*Ebeam + 27.0*(N - Z)/A + 0.4*Z*(A**(-1.0/3.0))
	vi = 0.0
	vsi = 13.5
	vso = 7.5
	vsoi = 0.0
	
	r0 = 1.25
	ri0 = 0.0
	rsi0 = 1.25
	rso0 = 1.25
	rsoi0 = 0.0
	
	a = 0.65
	ai = 0.0
	asi = 0.47
	aso = 0.47
	asoi = 0.0
	
	rc0 = 1.25
	
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Perey proton")
		PrintCalculatedQuantities(A,Z,E,Q)
	return string_list

# =============================================================================================== #
# Menet (protons)
def Menet(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
	
	# Calculate final parameters
	v = 49.9 - 0.22*E + 26.4*( N - Z )/A + 0.4*Z*( A**(-1.0/3.0) )
	vi = 1.2 + 0.09*E
	vsi = 4.2 - 0.05*E + 15.5*( N - Z )/A
	vso = 6.04
	vsoi = 0.0
	
	r0 = 1.16
	ri0 = 1.37
	rsi0 = 1.37
	rso0 = 1.064
	rsoi0 = 0.0
	
	a = 0.75
	ai = 0.74 - 0.008*E + float( N - Z )/float(A)
	asi = 0.74 - 0.008*E + float( N - Z )/float(A)
	aso = 0.78
	asoi = 00
	
	rc0 = 1.25
	
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Menet proton")
		PrintCalculatedQuantities(A,Z,E,Q)
	return string_list

 # =============================================================================================== #
# Varner (protons)
def Varner(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)

	# Calculate final parameters
	rc1 = 1.24*(A**(1.0/3.0) ) + 0.12
	ec = 1.73*Z/rc1
	eta =  float(N - Z)/float(A) 

	v = 52.9 + (13.1*float(N - Z)/float(A) ) + ( -0.299*( E - ec) )
	vi = 7.8/( 1 + np.exp( ( 35 - ( E - ec ) )/16.0 ) )
	vsi = ( 10 + ( 18.0*float(N - Z)/float(A) ) )/( 1 + np.exp( ( E - ec - 36.0 )/37.0 ) )
	vso = 5.9
	vsoi = 0.0
	
	r0 = ( ( 1.25*( A**(1.0/3.0) ) ) - 0.225 )*( A**(-1.0/3.0) )
	ri0 = ( ( 1.33*( A**(1.0/3.0) ) ) - 0.42 )*( A**(-1.0/3.0) )
	rsi0 = ( ( 1.33*( A**(1.0/3.0) ) ) - 0.42 )*( A**(-1.0/3.0) )
	rso0 = ( ( 1.34*( A**(1.0/3.0) ) ) - 1.2 )*( A**(-1.0/3.0) )
	rsoi0 = 0.0
	
	a = 0.69
	ai = 0.69
	asi = 0.69
	aso = 0.63
	asoi = 0.0

	rc0 = rc1*( A**(-1.0/3.0) )
	
	
	
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Varner proton")
		PrintCalculatedQuantities(A,Z,E,Q)
		print( "ec\t" + str(ec) )
		print( "eta\t" + str(eta) )
	return string_list

# =============================================================================================== #
# Becchetti and Greenlees (protons)
def BecchettiGreenlees(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# CHECK VALUE OF H
	CheckP(H)

	# CALCULATE TRIVIAL PARAMETERS
	[N, Q, E] = CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H)
	
	# Calculate final parameters
	v = 54.0 - 0.32*E + 0.4*Z*( A**(-1.0/3.0) ) + 24.0*( N - Z )/A
	vi = 0.22*E - 2.7 if 0.22*E - 2.7 > 0.0 else 0.0
	vsi = 11.8 - 0.25*E + 12.0*( N - Z )/A if 11.8 - 0.25*E + 12.0*( N - Z )/A > 0.0 else 0.0
	vso = 6.2
	vsoi = 0.0
	
	r0 = 1.17
	ri0 = 1.32
	rsi0 = 1.32
	rso0 = 1.01
	rsoi0 = 0.0
	
	a = 0.75
	ai = 0.51 + 0.7*( N - Z )/A
	asi = 0.51 + 0.7*( N - Z )/A
	aso = 0.75
	asoi = 0.0
	
	rc0 = 1.3
	
	# Format final paramaters into a list of strings
	v_list = [v, vi, vsi, vso, vsoi]
	r_list = [r0, ri0, rsi0, rso0, rsoi0]
	a_list = [a, ai, asi, aso, asoi]
	string_list = MakeStringList(v_list,r_list,a_list,rc0)
	
	if PRINT == 1:
		PrintOpticalModel(string_list, "Becchetti and Greenlees proton")
		PrintCalculatedQuantities(A,Z,E,Q)
	return string_list

# 124Te(p,d), Ebeam=15MeV, Ex=0MeV
#PRINT = 1
#BecchettiGreenlees(124, 52, 15, 0, 123.9028179, 1.00782503224, 2.01410177811, 122.9042698, 0)

