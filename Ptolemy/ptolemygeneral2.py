import pandas as pd
import os
import ptolemywriter as pt


#define a bunch of global variables which are parameters that will be fixed in the files
target = '116Cd'
A = float(target[0:3])
Z = 48
reaction = '(d,t)'
elab = 15

#give mass excess in MeV. get these off NNDC
delta_target = -88.7124765625
delta_projectile = 14.9498095703125
delta_ejectile = 13.1357216796875
delta_product = -88.0844765625

A_target = 116
A_projectile = 2
A_ejectile = 3
A_product = 115


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def masscalc(delta, A):
    return A + delta/931.5

M_Target = masscalc(delta_target, A_target)
M_Projectile = masscalc(delta_target, A_target)
M_Ejectile = masscalc(delta_target, A_target) 
M_Product = masscalc(delta_target, A_target)


#first set the location of the current directory, important as we want lots of changing directories
current_directory = os.getcwd()
savedir = os.getcwd() + '/'
ptolemydir = os.getcwd() + '/'
analysis_code_directory = current_directory[0:-7]

os.chdir('%sspectrum_analysis/output'%analysis_code_directory)
sa_dir = os.getcwd()

for root, dirs, filenames in os.walk(sa_dir):
    #filenames is a list/tuple of filenames
    print(filenames)
    peaks_df = pd.read_table(filenames[0], sep = ' ')

energylist = peaks_df['PEAK_ENERGY'].tolist()
energylist_mev = [0] * len(energylist)

for i in range(len(energylist)):
    energylist_mev[i] = energylist[i]/1000 
    
print('The energies to be calculated are:\n', energylist)

#go through states
for energy in energylist_mev:
    
    #get incoming potential
    deuteronomp = pt.AnCai(A, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)
    incoming_potential = ''
    for dparameter in deuteronomp:
        incoming_potential = incoming_potential + dparameter + '\n'

    
    #get outgoing potential
    protonomp = pt.koningDelaroche(A, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1)
    outgoing_potential = ''
    for pparameter in protonomp:
        outgoing_potential = outgoing_potential + pparameter + '\n'

    pt.ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir)


