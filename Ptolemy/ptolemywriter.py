import os
import sys
import numpy as np

# GLOBAL CONSTANTS
amu = 931.494	# amu in MeV/c^2
PRINT = 0
verbose = 0
sep = "\t"
div = "--------------------------------------------------"
DIV = "=================================================="

def koningDelaroche(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H1):
	# CHECK VALUE OF H1
	if H1 != 0 and H1 != 1:
		raise ValueError("p must have a value of 0 or 1")
	
	# CALCULATE TRIVIAL PARAMETERS
	# Number of neutrons
	N = A - Z
	
	# Calculate Q value
	Q = ( (M_Target + M_Projectile) - (M_Ejectile + M_Product) )*amu
	
	# Calculate energy
	if H1 == 0:
		E = Ebeam
	elif H1 == 1:
		E = Ebeam + Q - Ex
		
	if PRINT == 1:
		print(div)
		print("A" + sep + str(A))
		print("Z" + sep + str(Z))
		print("E" + sep + str(E))
		print("Q" + sep + str(Q))
	
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
	
	if PRINT == 1:
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
	stringList = []
	stringList.append("v = " + str(round(v,3)) + " r0 = " + str(round(r0,3)) + " a = " + str(round(a,3)))
	stringList.append("vi = " + str(round(vi,3)) + " ri0 = " + str(round(ri0,3)) + " ai = " + str(round(ai,3)))
	stringList.append("vsi = " + str(round(vsi,3)) + " rsi0 = " + str(round(rsi0,3)) + " asi = " + str(round(asi,3)))
	stringList.append("vso = " + str(round(vso,3)) + " rso0 = " + str(round(rso0,3)) + " aso = " + str(round(aso,3)))
	stringList.append("vsoi = " + str(round(vsoi,3)) + " rsoi0 = " + str(round(rsoi0,3)) + " asoi = " + str(round(asoi,3)) + " rc0 = " + str(round(rc0,3)))
	
	if verbose == 1:
		for i in range(0,len(stringList)):
			print(stringList[i])
	return stringList
	
def AnCai(A, Z, Ebeam, Ex, M_Target, M_Projectile, M_Ejectile, M_Product, H2):
	# CHECK VALUE OF H2 - Energy changes if ejectile is deuteron
	if H2 != 0 and H2 != 1:
		raise ValueError("p must have a value of 0 or 1")

	# CALCULATE TRIVIAL PARAMETERS
	# Number of neutrons
	N = A - Z
	
	# Calculate Q value
	Q = ( (M_Target + M_Projectile) - (M_Ejectile + M_Product) )*amu
	
	# Calculate energy
	if H2 == 0:
		E = Ebeam
	elif H2 == 1:
		E = Ebeam + Q - Ex
		
	if PRINT == 1:
		print(div)
		print("A" + sep + str(A))
		print("Z" + sep + str(Z))
		print("E" + sep + str(E))
		print(div)
		
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
	stringList = []
	stringList.append("v = " + str(round(v,3)) + " r0 = " + str(round(r0,3)) + " a = " + str(round(a,3)))
	stringList.append("vi = " + str(round(vi,3)) + " ri0 = " + str(round(ri0,3)) + " ai = " + str(round(ai,3)))
	stringList.append("vsi = " + str(round(vsi,3)) + " rsi0 = " + str(round(rsi0,3)) + " asi = " + str(round(asi,3)))
	stringList.append("vso = " + str(round(vso,3)) + " rso0 = " + str(round(rso0,3)) + " aso = " + str(round(aso,3)))
	stringList.append("vsoi = " + str(round(vsoi,3)) + " rsoi0 = " + str(round(rsoi0,3)) + " asoi = " + str(round(asoi,3)) + " rc0 = " + str(round(rc0,3)))
	
	if verbose == 1:
		for i in range(0,len(stringList)):
			print(stringList[i])

	return stringList


def ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir):


    #need reaction constituents, masses needed for calculating daughter nucleus
    #I'm keeping the 'H's there for now, because they might be necessary for optical model stuff I may or may not do later
    if reaction[1] == 'd':
        incomingparticle = '2H'
    elif reaction[1] == 'p':
        incomingparticle = '1H'

    if reaction[3] == 'd':
        outgoingparticle = '2H'
    elif reaction[3] == 'p':
        outgoingparticle = '1H'
    elif reaction[3] == 't':
        outgoingparticle = '3H'


    #get the daughter nucleus
    #pretty standard, chopping up strings to get a calculation for the mass, and since this is neutron transfer the element stays the same
    daughter = '%s%s'%(int(target[0:3])+int(incomingparticle[0])-int(outgoingparticle[0]), target[3:5])
    reaction_no_brackets = reaction[1:-1]
    e_kev = round(energy*1000,2)

    #make directories. One for all the files for this peak, one for input, output files
    directory_name = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
    os.mkdir('%s%s'%(savedir,directory_name))
    os.chdir('%s%s'%(savedir,directory_name))
    os.mkdir('%s%s/input_files'%(savedir,directory_name))
    os.mkdir('%s%s/output_files'%(savedir,directory_name))
    os.chdir(savedir)

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~WRITE PTOLEMY FILES AND RUN PTOLEMY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """

    #loop over l, l = 6 highest in shell model up to 126
    for l in range(0,7):

        #get parity, p = (-1)^l. This is important because If I don't have the parity right ptolemy will get annoyed. LAZY CODE can't even be bothered to work it out!
        p = None
        if l % 2 == 0: #even
            p = '+' 
        else: #odd
            p = '-'

        #need two mss as well! This is because ptolemy needs the j value
        for s in range(0,2):    

            #define j
            j = l + s - 0.5

            #don't want state of negative angular momentum!
            if j<0: continue

            #can't use this, needs to be in format x/2, also need to save a float j for later
            floatj = j
            j = '%s/2'%(int(j*2))

            #also needs to be for different principal quantum numbers. 
            for nodes in range(0,3):

                #pick the right n, l, j

                #first do this for p,d, neutron addition can populate the 50-82 or 82-126 subshell
                #check the nuclear shell model to see these states
                #we can't differentiate spins so I picked a random one for ls that could have either
                if incomingparticle == '1H' and outgoingparticle == '2H':

                    if l == 0 and (j != '1/2' or nodes != 2): continue                   # 2s1/2
                    if l == 1 and (j != '1/2' or nodes != 2): continue #j arbitrary      # 2p1/2, 2p3/2
                    if l == 2 and (j != '3/2' or nodes != 1): continue #j arbitrary      # 1d3/2, 1d5/2
                    if l == 3 and (j != '5/2' or nodes != 1): continue #j arbitrary      # 1f3/2, 1f5/2
                    if l == 4 and (j != '7/2' or nodes != 0): continue                   # 0g7/2
                    if l == 5 and (j != '11/2' or nodes != 0): continue #j arbitrary     # 0h9/2,0h11/2
                    if l == 6 and (j != '13/2' or nodes != 0): continue                  # 0i13/2 

                if (incomingparticle == '2H' and outgoingparticle == '1H') or (incomingparticle == '2H' and outgoingparticle == '3H'):
                    
                    if l == 0 and (j != '1/2' or nodes != 2): continue              # 2s1/2
                    if l == 1 and (j != '1/2' or nodes != 1): continue #j arbitrary # 1p1/2, 1p3/2
                    if l == 2 and (j != '3/2' or nodes != 1): continue #j arbitrary # 1d3/2, 1d5/2
                    if l == 3 and (j != '5/2' or nodes != 0): continue #j arbitrary # 0f5/2, 0f7/2
                    if l == 4 and (j != '7/2' or nodes != 0): continue #j arbitrary # 0g7/2, 0g9/2
                    if l == 5 and (j != '11/2' or nodes != 0): continue             # 0h11/2
                    if l == 6 : continue                                                # no l = 6 for (d,p)

                #write the ptolemy file using parameters defined so far
                ptolemyfile = """reset
r0target
print 0
REACTION: %s%s%s(%s%s %s) ELAB=%s
PARAMETERSET dpsb labangles r0target lstep=1 lmin=0 lmax=30 maxlextrap=0
LMIN=0 LMAX=40
PROJECTILE
NODES = 0
R = 1   A = 0.5   WAVEFUNCTION = av18   L = 0
ASYMPTOPIA = 500
;
TARGET
nodes=%s l=%s jp=%s r0=1.28 a=0.65 vso=6 rso0=1.10 aso=0.65 rc0=1.3 
;
INCOMING
%s;
OUTGOING
%s;
LABANGLES
ANGLEMIN=0 ANGLEMAX=60 ANGLESTEP=1
;
writens crosssec

end
"""%(target,reaction,daughter,j,p,energy,elab,nodes,l,j,incoming_potential,outgoing_potential)

                #now want to save this to a file in the new directory
                os.chdir('%s%s/input_files'%(savedir,directory_name))


                #need to define the name so we can write to different files for each parameter set
                infilename = 'Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s.in'%(target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p)

                #open the file and write to it
                f = open(infilename,'w')
                f.write(ptolemyfile)
                f.close()

                #ptolemy doesn't like brackets basically, so need to make backslashes in front of certain characters in some of the variables.

                directory_name2 = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
                infileptolemy = '<%s>'%infilename
                outfileptolemy = '%s%s/output_files/Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s.out'%(savedir,directory_name2,target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p)

                #run ptolemy for this file, specifying outfile path in the outfile directory
                os.system('%sptolemy %s %s'%(ptolemydir, infileptolemy, outfileptolemy))

                os.chdir(savedir)


    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PTCLEAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """



    #at this point all of the ptolemy output files are pretty messy, thankfully somebody's written a scipt to turn it into 2-column data
    #first of all go into the right directory, then run ptclean
    os.chdir('%s%s/output_files'%(savedir,directory_name))
    os.system('%sptclean'%ptolemydir)

    #want to remove the non-cleaned ones because they get in the way and are unnecessary. can remove this if we need them!
    #define the directory name, then loop over files in the directory
    indir = '%s%s/output_files'%(savedir,directory_name)
    for root, dirs, filenames in os.walk(indir):
        for f in filenames:
            #pick out the ones that don't say 'clean' at the end, and remove them
            if f[-5:-1] != 'clea': #yes I know it doesn't have the 'n'
                os.remove(f)
                

    return;







        
