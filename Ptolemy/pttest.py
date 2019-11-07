import numpy as np
import os, sys
import matplotlib.pyplot as plt
import glob

import ptolemywriter as pt
from opticalmodel_deuterons import *
from opticalmodel_protons import *

def potwriter(omp):
    pot = ''    
    for param in omp:
        pot = pot + param + '\n'
    return(pot)

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


current_directory = os.getcwd()
savedir = os.getcwd() + '/'
ptolemydir = os.getcwd() + '/'
analysis_code_directory = current_directory[0:-7]

os.chdir('%sspectrum_analysis/output'%analysis_code_directory)
sa_dir = os.getcwd()

files = glob.glob('*')

print(files)

xs = []
sxs = []
angles = []
energies = []

for fi in files:
    with open(fi) as f:
        header = f.readline()
        state1 = f.readline()
        state1 = state1.split()
        xs.append(float(state1[2]))
        sxs.append(float(state1[3]))
        angles.append(float(fi[16:18]))
        energies.append(float(state1[1]))


energy = np.average(np.array(energies))/1000
print(energy, energies)
xs = np.array(xs)
angles = np.array(angles)
sxs = np.array(sxs)

i = np.argsort(angles)

angles = angles[i]
xs = xs[i]
sxs = sxs[i]


os.chdir(ptolemydir)

inlist, outlist = [],[]

l = int(input("What is the angular momentum of this state"))

if reaction == '(d,p)':
    print('This is using potentials for a (d,p) reaction')
    #go through states

    inamelist = ["AnCai", "PereyPerey", "LohrHaeberli", "HanShiShen", "Bojowald", "DaehnickNR", "DaehnickR"]

    inlist.append(potwriter(AnCai(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(PereyPerey(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(LohrHaeberli(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(HanShiShen(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(Bojowald(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(DaehnickNR(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))
    inlist.append(potwriter(DaehnickR(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 0)))   


    onamelist = ["KoningDelaroche", "Perey", "Menet", "Varner", "BecchettiGreenlees"]

    outlist.append(potwriter(KoningDelaroche(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1))) 
    outlist.append(potwriter(Perey(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1))) 
    outlist.append(potwriter(Menet(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1))) 
    outlist.append(potwriter(Varner(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1))) 
    outlist.append(potwriter(BecchettiGreenlees(A_target, Z, elab, energy, M_Target, M_Projectile, M_Ejectile, M_Product, 1)))   


    for incoming_potential, inname in zip(inlist, inamelist):
        for outgoing_potential, outname in zip(outlist, onamelist):
            pt.pttestwrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir, l, inname, outname)

print(os.getcwd())

ofiles = glob.glob('*')
ofiles = np.sort(ofiles)

for file in ofiles: print(file)

sfs = []
fig = plt.figure(figsize = (11.69, 8.27))
labelax = fig.add_subplot(111)
labelax.spines['top'].set_color('none')
labelax.spines['bottom'].set_color('none')
labelax.spines['left'].set_color('none')
labelax.spines['right'].set_color('none')
labelax.tick_params(labelcolor = 'w', top = False, bottom = False, left = False, right = False)
labelax.set_xlabel('Angle (degrees)', size = 'xx-large', fontweight = 'bold')
labelax.set_ylabel('Cross Section (mb/sr)', size = 'xx-large', fontweight = 'bold')

nplots = len(sfs)
ncols = len(inamelist)
nrows = len(onamelist)

onamelist = ["Koning\nDelaroche", "Perey", "Menet", "Varner", "Becchetti\nGreenlees"]


for i, file in enumerate(ofiles):
    inp = file.split('_')[-3]
    outp = file.split('_')[-2]
    print(inp,outp)

    data = np.loadtxt(file)
    t_angles = data[:,0]
    t_xs = data[:,1]

    t_xsections_at_angles = np.array([np.zeros(len(angles))])
    numerator = 0
    denominator = 0

    #print(angles)
    #print(t_angles)

    for l in range(len(t_angles)):
        for m in range(len(angles)):
            #if there's a dodgy angle include this line
            #if angles[m] == 25: continue
            if angles[m] == t_angles[l] and xs[m] != 0: #don't want division by 0           
                numerator += xs[m]*t_xs[l] / sxs[m]**2
                denominator += t_xs[l]**2 / sxs[m]**2
                t_xsections_at_angles[0][m] = t_xs[l]
    try:
        norm = numerator/denominator
        norm_error = 1/np.sqrt(denominator)
    except:
        norm = 1
        print("warning. Normalisation set to 1 because denominator = 0")
        print(numerator, denominator)
    #turn the t_xsections into a numpy array to manipulate
    t_xsections_np = np.array(t_xs)
           
    #actually do the normalisation.
    #Divide the cross-section by the mean normalisation factor: sum_norm/number of angles
    norm_xsections = norm * t_xsections_np
    norm_xsections_at_angles = t_xsections_at_angles * norm

    sfs.append(xs[0]/t_xsections_at_angles[0][0])

    ax = fig.add_subplot(nrows, ncols, i+1)
    if i % ncols == ncols - 1:
        ax.set_ylabel(np.sort(onamelist)[int((i - ncols + 1)/ncols)])
        ax.yaxis.set_label_position('right')
    if i < 7:
        ax.set_xlabel(np.sort(inamelist)[i])    
        ax.xaxis.set_label_position('top') 
    ax.plot(t_angles, norm_xsections, color = 'xkcd:black')
    ax.errorbar(angles, xs, sxs, fmt = 'x', lw = 2)
    ax.set_xticks(np.arange(0,61,15))

fig.tight_layout()
plt.show()




sfs = np.array(sfs)
num_bins = 20
n, bins, patches = plt.hist(sfs, num_bins, range = (0,1), facecolor='blue', alpha=0.5)
plt.xlabel('Raw Spectroscopic Factor')
plt.ylabel('Frequency')
plt.show()

