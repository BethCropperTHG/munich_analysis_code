import pandas as pd
import os
import ptolemywriter as pt
from opticalmodel_deuterons import *
from opticalmodel_protons import  *
import glob


files = glob.glob('reaction_property_files/*')
print('Which parameter set would you like to use?' )
for i, f in enumerate(files):
    print(i, ":", f[-9:])

while True:
    i = input()
    if i.isdigit(): 
        i = int(i)
        break
    else: print('Warning, input invalid. Please try again.' )

with open(files[i]) as f:
    lines = f.readlines()
    target = lines[1][:-1]
    Z = int(lines[3][:-1])
    reaction = lines[5][:-1]
    elab = int(lines[7][:-1])
    
    delta_target = float(lines[10][:-1])
    delta_projectile = float(lines[12][:-1])
    delta_ejectile = float(lines[14][:-1])
    delta_product = float(lines[16][:-1])

    A_target = int(lines[18][:-1])
    A_projectile = int(lines[20][:-1])
    A_ejectile = int(lines[22][:-1])
    A_product = int(lines[24])

print("target = ", target)
print("Z = ", Z)
print("reaction = ", reaction)
print("elab = ", elab)
print("delta_target = " , delta_target)
print("delta_projectile = ", delta_projectile)
print("delta_ejectile = ", delta_ejectile)
print("delta_product = ", delta_product)
print("A_target = ", A_target)
print("A_projectile = ", A_projectile)
print("A_ejectile = ", A_ejectile)
print("A_product = ", A_product)



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

if reaction == '(d,p)':
    print('This is using potentials for a (d,p) reaction')
    #go through states
    for energy in energylist_mev:
    
        #get incoming potential
        deuteronomp = AnCai(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)
        incoming_potential = ''
        for dparameter in deuteronomp:
            incoming_potential = incoming_potential + dparameter + '\n'

    
        #get outgoing potential
        protonomp = KoningDelaroche(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1)
        outgoing_potential = ''
        for pparameter in protonomp:
            outgoing_potential = outgoing_potential + pparameter + '\n'

        pt.ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir)

if reaction == '(p,d)':
    print('This is using potentials for a (p,d) reaction')
    #go through states
    for energy in energylist_mev:
    
        #get outgoing potential
        deuteronomp = AnCai(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1)
        outgoing_potential = ''
        for dparameter in deuteronomp:
            outgoing_potential = outgoing_potential + dparameter + '\n'

    
        #get incoming potential
        protonomp = koningDelaroche(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)
        incoming_potential = ''
        for pparameter in protonomp:
            incoming_potential = incoming_potential + pparameter + '\n'
        
        pt.ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir)

if reaction == '(3HE,4HE)':
    print('This is using potentials for a (h,a) reaction')
    #go through states
    for energy in energylist_mev:
    
        #get outgoing potential
        alphaomp = BassaniPicard(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1)
        outgoing_potential = ''
        for aparameter in alphaomp:
            outgoing_potential = outgoing_potential + aparameter + '\n'

    
        #get incoming potential
        heliumomp = Pang(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)
        incoming_potential = ''
        for hparameter in heliumomp:
            incoming_potential = incoming_potential + hparameter + '\n'
        
        pt.ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir)
else:
    raise ValueError("You are trying to do calculations for a reaction you don't have the potentials for")




