#need to do a chi-squared analysis of each peak against predicted peaks from ptolemy, plotting the best assignments
#to match states to energies  order them (gs, 1st excited etc), then loop over the peaks in each spectrum and get chi-squared fits.plot the lowest one, or if the lowest two or 3 are close together, plot all of them and let the user know what's happened

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def spectroscopic_finder(exptdist, exptangles, t_dist, tangles, l, norm):
    
    #print('\n\n\n', exptdist, '\n\n\n', t_dist, '\n\n\n')

    '''
Experiment angles are 10, 18, 25, 31, 40 degrees. This is the peak of the l = 0,2,3,4,5
    '''
    spectroscopic_factor = None
    if l == 0:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 0)
    if l == 1:
        spectroscopic_factor = norm
    if l == 2:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 1)
    if l == 3:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 2)
    if l == 4:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 3)
    if l == 5:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 4)

    return spectroscopic_factor

def spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, angleindex):
    for i in range(len(tangles)):
        if exptangles[0] == tangles[i] and exptdist[0] != 0: #don't want division by 0
            spectroscopic_factor = exptdist[0]/t_dist[i] 

        elif exptangles[0] == tangles[i] and exptdist[0] == 0:
            spectroscopic_factor = norm
    
    return(spectroscopic_factor)

def intdiv(a, b):

    if a%b == 0:
        return a//b
    else:
        return a//b + 1

def colourplot(angle, peak_strength, peak_strengths_errors, t_dist, t_angle, l_selects, energy):
    
    final_dist = plt.gca()
    
    #final_dist.set_xlabel('Angle(degrees)')
    #final_dist.set_ylabel('Cross Section(mbarn)')


    if l_selects == '0':
        final_dist.plot(t_angle, t_dist, 'xkcd:black')
    elif l_selects == '1':
        final_dist.plot(t_angle, t_dist, 'xkcd:orange')
    elif l_selects == '2':
        final_dist.plot(t_angle, t_dist, 'xkcd:red')
    elif l_selects == '3':
        final_dist.plot(t_angle, t_dist, 'xkcd:brown')
    elif l_selects == '4':
        final_dist.plot(t_angle, t_dist, 'xkcd:green')
    elif l_selects == '5':  
        final_dist.plot(t_angle, t_dist, 'xkcd:blue')
    else:
        pass
    final_dist.errorbar(angle,peak_strength,peak_strengths_errors, fmt = 'x')
    #final_dist.legend(fontsize = 'small', numpoints = 0)
    plt.text(0.95,0.95, str(int(round(energy))) + ' keV', ha="right", va="top", transform=plt.gca().transAxes)

    return(final_dist)

def round_to_1(x):
   return round(x, -int(math.floor(math.log10(abs(x)))))

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


for i in range(len(state_dirs)):
    for j in range(len(state_dirs)-1-i):
        if float(state_dirs[j][83:]) > float(state_dirs[j+1][83:]):
            state_dirs[j], state_dirs[j+1] = state_dirs[j+1], state_dirs[j]
        pass

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET PEAKS AND ANGLES LISTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

plt.rc('axes', linewidth=3)
plt.rc('lines', linewidth=3)
plt.rc('xtick', labelsize = 'large')
plt.rc('xtick.major', width = 3)
plt.rc('ytick', labelsize = 'large')
plt.rc('ytick.major', width = 3)


os.chdir('%sspectrum_analysis/output'%analysis_code_directory)
sa_dir = os.getcwd()

spectra = []
angles = []
graphs = []

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
    l_list = []
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
            l = None
            
            if njp == "3_0.5+":
                l = 0
                handles.append("l = 0")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:black', lw = 3)
            if njp == "3_0.5-" or njp == "2_0.5-":
                l = 1
                handles.append("l = 1")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:orange', lw = 3)
            if njp == "2_1.5+":
                l = 2
                handles.append("l = 2")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:red', lw = 3)
            if njp == "1_2.5-" or njp == "2_2.5-":
                l = 3
                handles.append("l = 3")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:brown', lw = 3)
            if njp == "1_3.5+":
                l = 4
                handles.append("l = 4")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:green', lw = 3)
            if njp == "1_5.5-":
                l = 5
                handles.append("l = 5")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:blue', lw = 3)
            if njp == "1_6.5+":
                l = 6
                handles.append("l = 6")
                plt.plot(t_angles,norm_xsections[0], 'xkcd:purple', lw = 3)
            
            #plt.plot(t_angles, t_xsections)
     

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
            l_list.append(l)
            
            
            #print(handles)
            #print(peak_strengths)
            #print(norm_xsections)
            #print(chi_squared)

    index_min = min(range(len(chi_squared_list)), key=chi_squared_list.__getitem__)
    print('This peak has an energy of ', str(peak_energies[0]), 'keV\n') 
    print( 'The lowest chi-squared value is for a l = ' + str(l_list[index_min]), 'state\n')
    print( 'The chi-squared value for this l-assignment is ' + str(chi_squared_list[index_min]),)
    #legend(handles)
    #we want to keep other possible j-ps if they are close enough to the smallest one
    chi_squared_list_2 = [chi_squared_list, l_list, handles, n_dist_list, t_dist_list]
    
    #print('\n\n', chi_squared_list_2[0][1], '\n\n')
        
    for i in range(len(chi_squared_list_2[0])):
        for j in range(len(chi_squared_list)-1-i):
            if chi_squared_list_2[0][j] > chi_squared_list_2[0][j+1]:
                for lists in range(len(chi_squared_list_2)):
                    chi_squared_list_2[lists][j], chi_squared_list_2[lists][j+1] = chi_squared_list_2[lists][j+1], chi_squared_list_2[lists][j]
            pass
    

    are_there_other_ls = False
    has_the_message_been_displayed = False

    for chi2 in range(len(chi_squared_list)):

        if chi_squared_list_2[0][chi2] < 2 * chi_squared_list_2[0][0] and chi_squared_list_2[0][chi2] is not chi_squared_list_2[0][0]:
            if has_the_message_been_displayed == False:
                print('Other possible l assignments are:')
                has_the_message_been_displayed = True
            print('l = ', chi_squared_list_2[1][chi2], ' with a chi-squared of ', str(chi_squared_list_2[0][chi2]))

    print('\n\n\n')

    plt.errorbar(angles,peak_strengths,peak_strengths_error, fmt = 'x', lw = 2)
    
    
    fontsize = 14
    ax = plt.gca()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    #xlabel('Angle(Degrees)')
    #ylabel('Cross Section(mbarn)')
    plt.show()

#print(n_dist_list)
#steve = input('press any key to continue')
    
    

    #print(n_dist_list)
    #print(l_list)
    l_select = input('what l value does this state have?')    
    #l_select = '0'

    #lots of if statements to tell it which l to plot
    l_index = None
    for l in range(len(l_list)):
        if int(l_select) == l_list[l]:
            l_index = l

    colourplot(angles, peak_strengths, peak_strengths_error, n_dist_list[l_index][0], t_angles, l_select, peak_energies[0])

    spectroscopicFactor = spectroscopic_finder(peak_strengths, angles, t_dist_list[l_index], t_angles, int(l_select), t_dist_list[l_index][0]/n_dist_list[l_index][0])

    print('The spectroscopic factor for this state is:', spectroscopicFactor)
    
    #so you can't get the plots and simply paste them onto another set of axes, so we'll have to draw these axes again later. 
    dist_plotters = [angles, peak_strengths, peak_strengths_error, n_dist_list[l_index][0], t_angles, l_select, peak_energies[0]]
    graphs.append(dist_plotters)
    plt.show()  

#print(graphs[0])
#print(peak_energies)

#get how many rows and columns there are
nplots = len(graphs)

npages = intdiv(nplots, 24)

for pages in range(npages):

    ncols = 4
    nrows = 6

    #make an a4 figure
    fig = plt.figure(figsize = (8.27, 11.69))
    
    labelax = fig.add_subplot(111)

    #loop over the graphs and add the subplots
    for i in range(24):
        ax = fig.add_subplot(nrows, ncols, i+1)
        #print('\n\n\n', graphs[i][4], '\n\n\n', graphs[i][3])
        try:
            ax = colourplot(graphs[i + pages * 24][0], graphs[i + pages * 24][1],graphs[i + pages * 24][2],graphs[i + pages * 24][3], graphs[i + pages * 24][4], graphs[i + pages * 24][5],graphs[i + pages * 24][6] )
        except:
            ax.axis('off')
            break
        plt.xticks(np.arange(0,61,15))
    
        if max(graphs[i + pages * 24][1]) > max(graphs[i + pages * 24][3]):
            plt.yticks(np.arange(0, max(graphs[i + pages * 24][1])+max(graphs[i + pages * 24][1])/10, round_to_1(max(graphs[i + pages * 24][1])/5)))
        else:
            plt.yticks(np.arange(0, max(graphs[i + pages * 24][3])+max(graphs[i + pages * 24][3])/10, round_to_1(max(graphs[i + pages * 24][3])/5)))

    
        if i < 20:
            plt.setp(ax.get_xticklabels(), visible = False)

    '''    
    if i == 21:
        ax.set_xlabel('Angle(degrees)', size = 'large', weight = 'bold')

    if i == 8:
        ax.set_ylabel('Cross-section(mb/sr)', size = 'large', weight = 'bold')
    '''        


    fig.align_ylabels(axs=None)

    
    labelax.spines['top'].set_color('none')
    labelax.spines['bottom'].set_color('none')
    labelax.spines['left'].set_color('none')
    labelax.spines['right'].set_color('none')
    labelax.tick_params(labelcolor = 'w', top = False, bottom = False, left = False, right = False)
    labelax.set_xlabel('Angle (degrees)', size = 'large', fontweight = 'bold')
    labelax.set_ylabel('Cross Section (mb/sr)', size = 'large', fontweight = 'bold')

    fig.tight_layout()
    plt.show()

