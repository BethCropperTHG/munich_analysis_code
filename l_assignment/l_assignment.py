#need to do a chi-squared analysis of each peak against predicted peaks from ptolemy, plotting the best assignments
#to match states to energies  order them (gs, 1st excited etc), then loop over the peaks in each spectrum and get chi-squared fits.plot the lowest one, or if the lowest two or 3 are close together, plot all of them and let the user know what's happened

import os
import pandas as pd
import numpy as np
from matplotlib.pyplot import *


#first set the location of the current directory, important as we want lots of changing directories
current_directory = os.getcwd()
analysis_code_directory = current_directory[0:-12]

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET PTOLEMY OUTFILE DIRECTORY LIST THEN SORT IT~~~~~~~~~~~~~~~~~~~~
'''

#first thing to do is get a list of directories that contain ptolemy outfiles, in ascending order of energy
#change to the ptolemy directory, and create an empty list which will contain the directories
os.chdir('%s/Ptolemy'%analysis_code_directory)
ptolemydir = os.getcwd()
state_dirs = []

#loop over the directories in the ptolemy folder and root out the ones that contain an output_files subdirectory
for dirs, subdirs, filenames in os.walk(ptolemydir):
    #need a try block here since some directories won't have subdirectories. we want to ignore those so do nothing if that 
    #exception has been encountered
 
    try:
        if 'output_files' in subdirs:
            state_dirs.append(dirs)
    except:
        pass          
#now need to sort this list. need to put the numbers in ascending order. to do this i can't use the built in sort.
#this is because i'm sorting the numbers in a big string. this means I'll have to write my own sorting algorithm
#bubble sort is easy and shouldn't be too slow for this purpose
for d in state_dirs:
    print(d)

for i in range(len(state_dirs)):
    for j in range(len(state_dirs)-1-i):
        if float(state_dirs[j][83:]) > float(state_dirs[j+1][83:]):
            state_dirs[j], state_dirs[j+1] = state_dirs[j+1], state_dirs[j]
        pass

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET PEAKS AND ANGLES LISTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

rc('axes', linewidth=2)
os.chdir('%sspectrum_analysis/output'%analysis_code_directory)
sa_dir = os.getcwd()

spectra = []
angles = []

for root, dirs, filenames in os.walk(sa_dir):
    #f is a string of a filename, and filenames is a list/tuple of filenames
    for f in filenames:
        angle = float(f[16:-24])
        peaks_df = pd.read_table(f, sep = ' ')
        #print(peaks_df)
        spectra.append(peaks_df)
        angles.append(angle)
        
        
peaks_no = peaks_df.shape[0]

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET THEORETICAL AND EXPERIMENTAL ANGULAR DISTRIBUTIONS~~~~~~~~~~~~
'''

#loop over different peaks, want to do this for all of them!
for i in range(peaks_no):
    #first get experimental angular distribution
    #make empty lists for expt angular distribution and uncertainty
    peak_strengths = []
    peak_energies = []
    peak_strengths_error = []
    #loop over the spectra to find the cross section for peak i at different angles
    
    for j in range(len(angles)):
        #find the cross sections and errors and stick them on the angular distribution list
        xsection = spectra[j].CROSS_SECTION[i]
        peak_strengths.append(xsection)

        sxsection = spectra[j].ERROR[i]
        peak_strengths_error.append(sxsection)

        energy = spectra[j].PEAK_ENERGY[i]
        peak_energies.append(energy)     
        
    #now we need to get the correspoding theoretical angular distribution.
    #THIS USES THE SORTED LIST FROM EARLIER TO GET THE THEORETICAL DISTRIBUTION
    theor_dist_dir = '%s/output_files'%state_dirs[i]
    os.chdir(theor_dist_dir)
    #print('\n\n\n')
    #print(theor_dist_dir)
    #print('\n\n\n')

    #now we need to loop over the theoretical distributions and see which ones fit best
    #make a list for the chi-squareds
    chi_squared_list = []
    jp_list = []
    handles = []
    t_dist_list = []
    n_dist_list = []    
    
    #Get the theoretical distributions
    for root, dirs, filenames in os.walk(theor_dist_dir):
        for file in filenames:
            #Make empty lists for the angles and cross sections to be stuck on from the file
            t_angles = []
            t_xsections = []
            f = open(file)
            for line in f:
                splitline = line.split()
                t_angles.append(float(splitline[0]))
                t_xsections.append(float(splitline[1]))
            f.close()
            
            #Now we do a normalisation of the theoretical data to the experimental one
            
            #make an empty list to stick the normalised cross sections on.
            #In this case we use a numpy array since we want to do maths to get this one
            t_xsections_at_angles = np.array([np.zeros(len(angles))])
            #also want to initialise the normalisation factor
            numerator = 0
            denominator = 0

            #loop over the theoretical angles and experimental ones to find the matching ones
            #average ratio of experimental values and theoretical ones is the normlisation
            #for each matching value, this loop calculates the normalisation factor and adds it to the sum_norm value
            for l in range(len(t_angles)):
                for m in range(len(angles)):
                    if angles[m] == t_angles[l] and peak_strengths[m] != 0: #don't want division by 0
                        numerator += peak_strengths[m]*t_xsections[l] / peak_strengths_error[m]**2
                        denominator += t_xsections[l]**2 / peak_strengths_error[m]**2
                        t_xsections_at_angles[0][m] = t_xsections[l]
            try:
                norm = numerator/denominator
            except:
                norm = 1
            #turn the t_xsections into a numpy array to manipulate
            t_xsections_np = np.array([t_xsections])
            
            #actually do the normalisation.
            #Divide the cross-section by the mean normalisation factor: sum_norm/number of angles
            norm_xsections = norm * t_xsections_np
            norm_xsections_at_angles = t_xsections_at_angles * norm



            #plot this so I can visualise the normalised fits against the experimental data

            njp = file[-24:-23] + '_' + file[-14:-10]
            
            if njp == "3_0.5+":
                handles.append("l = 0")
                plot(t_angles,norm_xsections[0], 'xkcd:black', lw = 3)
            if njp == "3_0.5-" or njp == "2_0.5-":
                handles.append("l = 1")
                plot(t_angles,norm_xsections[0], 'xkcd:orange', lw = 3)
            if njp == "2_1.5+":
                handles.append("l = 2")
                plot(t_angles,norm_xsections[0], 'xkcd:red', lw = 3)
            if njp == "1_2.5-" or njp == "2_2.5-":
                handles.append("l = 3")
                plot(t_angles,norm_xsections[0], 'xkcd:brown', lw = 3)
            if njp == "1_3.5+":
                handles.append("l = 4")
                plot(t_angles,norm_xsections[0], 'xkcd:green', lw = 3)
            if njp == "1_5.5-":
                handles.append("l = 5")
                plot(t_angles,norm_xsections[0], 'xkcd:blue', lw = 3)
            if njp == "1_6.5+":
                handles.append("l = 6")
                plot(t_angles,norm_xsections[0], 'xkcd:purple', lw = 3)
            
            #plot(t_angles, t_xsections)
     

            #I also need to save this so I can do a single plot
            t_dist_list.append(t_xsections)
            n_dist_list.append(norm_xsections)

            #Time for chi-squared analysis

            chi_squared = 0

            #chi-squared = S((O-E)^2/sigma^2)
            for m in range(len(angles)):
                if peak_strengths[m] != 0:
                    #work out the contribution to the chi squared
                    residual = peak_strengths[m] - norm_xsections_at_angles[0][m]
                    chi_squared_contribution = residual**2/peak_strengths_error[m]**2
                    chi_squared += chi_squared_contribution

            
            
            
            #handles.append(njp)
            chi_squared_list.append(chi_squared)
            jp_list.append(njp)
            
            #print(handles)
            #print(peak_strengths)
            #print(norm_xsections)
            #print(chi_squared)

    index_min = min(range(len(chi_squared_list)), key=chi_squared_list.__getitem__)
    print('This peak has an energy of ' + str(peak_energies[0])) 
    print( 'keV. The lowest chi-squared value is for a ' + str(jp_list[index_min]))
    print( 'state. The chi-squared value for this l-assignment is ' + str(chi_squared_list[index_min]))
    #legend(handles)
    #we want to keep other possible j-ps if they are close enough to the smallest one
    chi_squared_list_2 = [chi_squared_list, jp_list, handles]
    
        
    for i in range(len(chi_squared_list_2[0])):
        for j in range(len(chi_squared_list)-1-i):
            if chi_squared_list_2[0][j] > chi_squared_list_2[0][j+1]:
                chi_squared_list_2[0][j], chi_squared_list_2[0][j+1] = chi_squared_list_2[0][j+1], chi_squared_list_2[0][j]
                chi_squared_list_2[1][j], chi_squared_list_2[1][j+1] = chi_squared_list_2[1][j+1], chi_squared_list_2[1][j]
                chi_squared_list_2[2][j], chi_squared_list_2[2][j+1] = chi_squared_list_2[2][j+1], chi_squared_list_2[2][j]
            pass
    

    are_there_other_ls = False
    has_the_message_been_displayed = False

    for chi2 in range(len(chi_squared_list)):

        if chi_squared_list_2[0][chi2] < 2 * chi_squared_list_2[0][0] and chi_squared_list_2[0][chi2] is not chi_squared_list_2[0][0]:
            if has_the_message_been_displayed == False:
                print('Other possible l assignments are:')
                has_the_message_been_displayed = True
            print('njp = ' + chi_squared_list_2[1][chi2] +' with a chi-squared of ' + str(chi_squared_list_2[0][chi2]))

    errorbar(angles,peak_strengths,peak_strengths_error, fmt = 'x', lw = 2)
    
    
    fontsize = 14
    ax = gca()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    #xlabel('Angle(Degrees)')
    #ylabel('Cross Section(mbarn)')
    show()

#print(n_dist_list)
#steve = input('press any key to continue')
    
#    final_dist = figure()
#    final_dist = final_dist.add_subplot(111)
#    final_dist.set_xlabel('Angle(degrees)')
#    final_dist.set_ylabel('Cross Section(mbarn)')

    #print(n_dist_list)
#    l_select = input('what l value does this state have?')    

    #lots of if statements to tell it which l to plot
#    if l_select == '0':
#        final_dist.plot(t_angles, n_dist_list[6][0], 'xkcd:black')
#    elif l_select == '1':
#        final_dist.plot(t_angles, n_dist_list[5][0], 'xkcd:orange')
#    elif l_select == '2':
#        final_dist.plot(t_angles, n_dist_list[3][0], 'xkcd:red')
#    elif l_select == '3':
#        final_dist.plot(t_angles, n_dist_list[4][0], 'xkcd:brown')
#    elif l_select == '4':
#        final_dist.plot(t_angles, n_dist_list[0][0], 'xkcd:green')
#    elif l_select == '5':  
#        final_dist.plot(t_angles, n_dist_list[1][0], 'xkcd:blue')
#    else:
#        pass
#    final_dist.errorbar(angles,peak_strengths,peak_strengths_error, fmt = 'x')
    

#    show()  


