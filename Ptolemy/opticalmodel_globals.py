# opticalmodel_globals
# Contains global variables for the optical model calculations
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# =============================================================================================== #
# GLOBAL CONSTANTS
amu = 931.494	# amu in MeV/c^2
PRINT = 0
sep = "\t"
div = "--------------------------------------------------"
DIV = "=================================================="
# ----------------------------------------------------------------------------------------------- #
# Print out the calculated quantities in the string list
def PrintOpticalModel(string_list, name):
	print(DIV)
	print(name + " potential:")
	print(div)
	for i in range(0,len(string_list)):
			print(string_list[i])

# ----------------------------------------------------------------------------------------------- #
def PrintCalculatedQuantities(A,Z,E,Q):
	print(div)
	print("A" + sep + str(A))
	print("Z" + sep + str(Z))
	print("E" + sep + str(E))
	print("Q" + sep + str(Q))

# ----------------------------------------------------------------------------------------------- #
# Calculate Q value
def CalcQ( M_Target, M_Projectile, M_Ejectile, M_Product ):
	return ( (M_Target + M_Projectile) - (M_Ejectile + M_Product) )*amu

# ----------------------------------------------------------------------------------------------- #
# Calculate trivial quantities
def CalcTrivials(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H):
	# Number of neutrons
	N = A - Z
	
	# Calculate Q value
	Q = CalcQ( M_Target, M_Projectile, M_Ejectile, M_Product )
	
	# Calculate energy
	if H == 0:
		E = Ebeam
	elif H == 1:
		E = Ebeam + Q - Ex

	# Return list
	triv_list = [N, Q, E]

	return triv_list

# ----------------------------------------------------------------------------------------------- #
# Make the string list from the calculated quantities
def MakeStringList(v,r,a,rc0):
	stringList = []
	stringList.append("v = " +    str(round(v[0],3)) + " r0 = " +    str(round(r[0],3)) + " a = " +    str(round(a[0],3)))
	stringList.append("vi = " +   str(round(v[1],3)) + " ri0 = " +   str(round(r[1],3)) + " ai = " +   str(round(a[1],3)))
	stringList.append("vsi = " +  str(round(v[2],3)) + " rsi0 = " +  str(round(r[2],3)) + " asi = " +  str(round(a[2],3)))
	stringList.append("vso = " +  str(round(v[3],3)) + " rso0 = " +  str(round(r[3],3)) + " aso = " +  str(round(a[3],3)))
	stringList.append("vsoi = " + str(round(v[4],3)) + " rsoi0 = " + str(round(r[4],3)) + " asoi = " + str(round(a[4],3)) + " rc0 = " + str(round(rc0,3)))
	return stringList

# ----------------------------------------------------------------------------------------------- #
def CheckP(p):
	if p != 0 and p != 1:
		raise ValueError("p must have a value of 0 or 1")







































