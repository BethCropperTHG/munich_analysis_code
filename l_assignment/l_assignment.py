#need to do a chi-squared analysis of each peak against predicted peaks from ptolemy, plotting the best assignments
#to match states to energies  order them (gs, 1st excited etc), then loop over the peaks in each spectrum and get chi-squared fits.plot the lowest one, or if the lowest two or 3 are close together, plot all of them and let the user know what's happened

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def ticklengthsetter(dist, errors = None, ticksno = 2):
    if errors is None:
        maximum = max(dist)
    else:
        maximum = max(dist) + max(errors)

    magnitude = math.floor(math.log10(maximum))
    firstnumber = first_number(maximum)
    maxtick = (firstnumber + 1 ) * 10 ** magnitude
    tickstep = maxtick/ticksno
    return(np.arange(0, maxtispectroscopic_finderck + 0.000001, tickstep))

def first_number(number):
    numberstring = str(number)
    for digits in numberstring:
        if digits is not '0':
            if digits is not '.':
                return float(digits)
                break 

def spectroscopic_finder(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm):
    
    #print('\n\n\n', exptdist, '\n\n\n', t_dist, '\n\n\n')

    
#Experiment angles are 10, 18, 25, 31, 40 degrees. for (d,p) This is the peak of the l = 0,2,3,4,5
#For (p,d), it's at 8, 17, 26, 31, 39, l = 0, 2, 3, 4, 5
    '''
    spectroscopic_factor = None
    if l == 0:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 0)
        #spectroscopic_factor = norm
    if l == 1:
        spectroscopic_factor = norm
        #spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 0)
    if l == 2:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 1)
    if l == 3:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 2)
    if l == 4:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 3)
    if l == 5:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 4)
    if l == 6:
        spectroscopic_factor = norm

    return spectroscopic_factor
    '''
    '''
    spectroscopic_factor = None
    if l == 0:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 0)
        #spectroscopic_factor = norm
    if l == 1:
        spectroscopic_factor = norm
        #spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 1)
    if l == 2:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 1)
    if l == 3:
        #spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 2)
        spectroscopic_factor = norm
    if l == 4:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 3)
    if l == 5:
        spectroscopic_factor = spectroscopic_calculator(exptdist, exptangles, t_dist, tangles, l, norm, 4)
    '''
    
    spectroscopic_factor = [None, None]
    if l == 0:
        spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 0)
        #spectroscopic_factor[0] = norm
        #spectroscopic_factor[1] = err_norm
    if l == 1:
        spectroscopic_factor[0] = norm
        spectroscopic_factor[1] = err_norm
        #spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 1)
    if l == 2:
        spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 1)
    if l == 3:
        spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 2)
        #spectroscopic_factor[0] = norm
        #spectroscopic_factor[1] = err_norm
    if l == 4:
        spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 3)
    if l == 5:
        spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 4)
    
    #spectroscopic_factor = [None, None]
    #spectroscopic_factor = spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, 0)


    return spectroscopic_factor

def spectroscopic_calculator(exptdist, sexptdist, exptangles, t_dist, tangles, l, norm, err_norm, angleindex):
    spectroscopic_factor = [None,None]

    sortedangles = np.sort(exptangles)
    used_angle = sortedangles[angleindex]

    for i in range(len(tangles)):
        if used_angle == tangles[i] and exptdist[int(np.where(exptangles == used_angle)[0])] != 0: #don't want division by 0
            #print(exptdist, exptangles, t_dist[i])
            #plt.plot(exptangles,exptdist, 'o')
            #plt.plot(tangles,t_dist)
            #plt.show()
            spectroscopic_factor[0] = exptdist[int(np.where(exptangles == used_angle)[0])]/t_dist[i]
            spectroscopic_factor[1] = sexptdist[int(np.where(exptangles == used_angle)[0])]/t_dist[i]

        elif used_angle == tangles[i] and exptdist[int(np.where(exptangles == used_angle)[0])] == 0:
            spectroscopic_factor[0] = norm
            spectroscopic_factor[1] = err_norm
            print('Warning: This spectroscopic factor has been obtained from the normalisation rather than the peak')
    
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
        final_dist.plot(t_angle, t_dist, 'xkcd:light brown')
    elif l_selects == '4':
        final_dist.plot(t_angle, t_dist, 'xkcd:green')
    elif l_selects == '5':  
        final_dist.plot(t_angle, t_dist, 'xkcd:blue')
    elif l_selects == '6':  
        final_dist.plot(t_angle, t_dist, 'xkcd:purple')
    else:
        return None
    
    #create angle distribution which doesn't plot the zeroes
    new_angle = []
    new_peak_strength = []
    new_peak_strength_error = []
    for i in range(len(peak_strength)):
        if peak_strength[i] != 0:
            new_angle.append(angle[i])
            new_peak_strength.append(peak_strength[i])
            new_peak_strength_error.append(peak_strengths_errors[i])

    final_dist.errorbar(new_angle,new_peak_strength,new_peak_strength_error, color = 'xkcd:black', fmt = '.')
    #final_dist.legend(fontsize = 'small', numpoints = 0)
    plt.text(0.95,0.95, str(int(round(energy))) + ' keV', ha="right", va="top", transform=plt.gca().transAxes)

    return(final_dist)

def colourplot_doublet(angle, peak_strength, peak_strength_error, dist1, l1, dist2, l2, t_angle, energy):

    ovr_fig = colourplot(angle, peak_strength, peak_strength_error, dist1, t_angle, l1, energy)

    if l2 == '0':
        ovr_fig.plot(t_angle, dist2, 'xkcd:black')
    elif l2 == '1':
        ovr_fig.plot(t_angle, dist2, 'xkcd:orange')
    elif l2 == '2':
        ovr_fig.plot(t_angle, dist2, 'xkcd:red')
    elif l2 == '3':
        ovr_fig.plot(t_angle, dist2, 'xkcd:light brown')
    elif l2 == '4':
        ovr_fig.plot(t_angle, dist2, 'xkcd:green')
    elif l2 == '5':  
        ovr_fig.plot(t_angle, dist2, 'xkcd:blue')
    elif l2 == '6':  
        ovr_fig.plot(t_angle, dist2, 'xkcd:purple')
    else:
        pass    

    ovr_fig.plot(t_angles, dist1 + dist2, ls = 'dashed', color = 'xkcd:dark teal')

    return(ovr_fig)


def round_to_1(x):
   return round(x, -int(math.floor(math.log10(abs(x)))))

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i
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


'''
#if autofit
for i in range(len(state_dirs)):
    for j in range(len(state_dirs)-1-i):
        if float(state_dirs[j][87:]) < float(state_dirs[j+1][87:]):
            state_dirs[j], state_dirs[j+1] = state_dirs[j+1], state_dirs[j]
        pass
'''

#if manual
for i in range(len(state_dirs)):
    for j in range(len(state_dirs)-1-i):
        if float(state_dirs[j][83:]) > float(state_dirs[j+1][83:]):
            state_dirs[j], state_dirs[j+1] = state_dirs[j+1], state_dirs[j]
pass

print(len(state_dirs))

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET PEAKS AND ANGLES LISTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

plt.rc('axes', linewidth=3)
plt.rc('lines', linewidth=3)
plt.rc('xtick', labelsize = 'large')
plt.rc('xtick.major', width = 3)
plt.rc('xtick.minor', width = 2, size = 3)
plt.rc('ytick', labelsize = 'large')
plt.rc('ytick.major', width = 3)
plt.rc('ytick.minor', width = 2, size = 3)


os.chdir('%sspectrum_analysis/output'%analysis_code_directory)
sa_dir = os.getcwd()
 

spectra = []
angles = []
graphs = []
reactionname = None

#create pandas df for spectroscopic factors to be added to
spectroscopic_df = pd.DataFrame(data = None, columns = ['ENERGY', 'L', 'SPECTROSCOPIC_FACTOR', 'ERROR'])

for root, dirs, filenames in os.walk(sa_dir):
    #f is a string of a filename, and filenames is a list/tuple of filenames
    print('\n\n', filenames, '\n\n')
    for f in filenames:
        reactionname = f[0:16] + f[26:32]
        print('The reaction is: ', reactionname)
        print(f) 
        angle = float(f[16:-24])
        print(angle, '\n')
        peaks_df = pd.read_table(f, sep = ' ')
        #print(peaks_df)
        spectra.append(peaks_df)
        angles.append(angle)
        
        
peaks_no = peaks_df.shape[0]
print(peaks_no)
'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET THEORETICAL AND EXPERIMENTAL ANGULAR DISTRIBUTIONS~~~~~~~~~~~~
'''

#loop over different peaks, want to do this for all of them!
for i in range(peaks_no):

    #first get experimental angular distribution
    #make empty lists for expt angular distribution and uncertainty
    peak_strengths = []
    peak_energies = []
    peak_positions = []
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
         
        position = spectra[j].PEAK_POSITION[i]
        peak_positions.append(position)      
    peak_energy = np.mean(np.array(peak_energies))    
    #now we need to get the correspoding theoretical angular distribution.
    #THIS USES THE SORTED LIST FROM EARLIER TO GET THE THEORETICAL DISTRIBUTION

    print(i)
    print(state_dirs[i])    
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
    n_dist_list_at_angles = []
    norm_errors_list = [] 

    #Get the theoretical distributions
    for root, dirs, filenames in os.walk(theor_dist_dir):
        for file in filenames:
            if (file[-24:-23] + '_' + file[-14:-10]) == "2_2.5+": continue
            #Make empty lists for the angles and cross sections to be stuck on from the file
            t_angles = []
            t_xsections = []
            f = open(file)
            #print('The theoretical distribution file is: ',f)
            for i, line in enumerate(f):
                #print('i = ', i, 'length = ', file_len(file))
                if i is not file_len(file):
                    splitline = line.split()
                    t_angles.append(float(splitline[0]))
                    t_xsections.append(float(splitline[1]))
            f.close()
            #print(len(t_angles))
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
                    #if there's a dodgy angle include this line
                    #if angles[m] == 25: continue
                    if angles[m] == t_angles[l] and peak_strengths[m] != 0: #don't want division by 0           
                        numerator += peak_strengths[m]*t_xsections[l] / peak_strengths_error[m]**2
                        denominator += t_xsections[l]**2 / peak_strengths_error[m]**2
                        t_xsections_at_angles[0][m] = t_xsections[l]
            try:
                norm = numerator/denominator
                norm_error = 1/math.sqrt(denominator)
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
                plt.plot(t_angles,norm_xsections[0], 'xkcd:light brown', lw = 3)
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
            n_dist_list_at_angles.append(norm_xsections_at_angles)
            norm_errors_list.append(norm_error)

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
    chi_squared_list_2 = [chi_squared_list, l_list, handles, n_dist_list, t_dist_list, n_dist_list_at_angles, norm_errors_list]
    
    #print('\n\n', chi_squared_list_2[0][1], '\n\n')
    '''    
    for i in range(len(chi_squared_list_2[0])):
        for j in range(len(chi_squared_list)-1-i):
            if chi_squared_list_2[0][j] > chi_squared_list_2[0][j+1]:
                for lists in range(len(chi_squared_list_2)):
                    chi_squared_list_2[lists][j], chi_squared_list_2[lists][j+1] = chi_squared_list_2[lists][j+1], chi_squared_list_2[lists][j]
            pass
    '''
    for i in range(len(chi_squared_list_2[0])):
        for j in range(len(chi_squared_list)-1-i):
            if chi_squared_list_2[0][j] > chi_squared_list_2[0][j+1]:
                for lists in range(len(chi_squared_list_2)):
                    #print('j = ', j, 'max_j = ',range(len(chi_squared_list)-1-i), 'lists = ', lists, 'max_lists = ',len(chi_squared_list_2))
                    #print('Things we are trying to assign are:')
                    #print(chi_squared_list_2[lists][j])
                    #print(chi_squared_list_2[lists][j+1])                    
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
    #print(l_list)
    l_select = input('What l value does this state have? Input \'d\' if the state is a doublet, and \'n\' if you don\'t want to assign it')    
    #l_select = '0'


    l_index = None

    try:
        for l in range(len(l_list)):
            if int(l_select) == l_list[l]:
                l_index = l
    except ValueError:
        pass

    #now to deal with doublets
    if l_select == 'd':
        ls = input('Input the ls that this state has, with no spaces')
        l_1 = ls[0]
        l_2 = ls[1]

        l1_index = None
        l2_index = None

        #get the indices for these ls        
        for l in range(len(l_list)):
            if int(l_1) == l_list[l]:
                l1_index = l

        for l in range(len(l_list)):
            if int(l_2) == l_list[l]:
                l2_index = l
        
        #Work out weights for two theoretical distributions
        dist1 = np.array(n_dist_list_at_angles[l1_index][0])
        dist2 = np.array(n_dist_list_at_angles[l2_index][0])

        full_dist1 = np.array(n_dist_list[l1_index][0])
        full_dist2 = np.array(n_dist_list[l2_index][0])
        
        target = np.array(peak_strengths)
        s_target = np.array(peak_strengths_error)

        npangles = np.array(angles)

        s_target = s_target[np.where(target != 0)]
        npangles = npangles[np.where(target != 0)]
        target = target[np.where(target != 0)]


        func1 = np.polyfit(angles, dist1, len(target) - 1)
        func2 = np.polyfit(angles, dist2, len(target) - 1)
        

        #print(target)
        #print(s_target)
        #print(dist1)

        weights_guess = 0.5
        bnds = (0.,1.)

        def objective(x, A):
            if len(target) == 5:
                return A * (func1[0] * x**4 + func1[1] * x**3 + func1[2] * x**2 + func1[3] * x + func1[4]) + (1-A) * (func2[0] * x**4 + func2[1] * x**3 + func2[2] * x**2 + func2[3] * x + func2[4])
            if len(target) == 4:
                return A * (func1[0] * x**3 + func1[1] * x**2 + func1[2] * x**1 + func1[3] * x**0) + (1 - A) * (func2[0] * x**3 + func2[1] * x**2 + func2[2] * x**1 + func2[3] * x**0)
        
        '''
        x = np.array(t_angles)        
        #print(t_angles, func1[0] * x**4 + func1[1] * x**3 + func1[2] * x**2 + func1[3] * x + func1[4]) + (1-A) * (func2[0] * x**4 + func2[1] * x**3 + func2[2] * x**2 + func2[3] * x + func2[4])
        plt.plot(t_angles, func1[0] * x**4 + func1[1] * x**3 + func1[2] * x**2 + func1[3] * x + func1[4])
        plt.plot(t_angles, func2[0] * x**4 + func2[1] * x**3 + func2[2] * x**2 + func2[3] * x + func2[4])
        plt.plot(angles, dist1, 'o')
        plt.plot(angles, dist2, 'o')
        #plt.errorbar(angles, target, s_target, fmt = '.')
        plt.show()
        '''

        import scipy.optimize as optimization
        optimised = optimization.curve_fit(objective, npangles, target, weights_guess, s_target, bounds = bnds)

        sumdist = optimised[0][0] * full_dist1 + (1 - optimised[0][0]) * full_dist2

        print(optimised)
        '''        
        plt.plot(t_angles, sumdist)
        plt.plot(t_angles, optimised[0][0] * full_dist1)
        plt.plot(t_angles, optimised[0][1] * full_dist2)        
        plt.errorbar(angles,target,s_target)
        plt.show()
        '''


        #do a colourplot with a second theoretical distribution in there
        colourplot(angles, peak_strengths, peak_strengths_error, optimised[0][0] * full_dist1, t_angles, l_1, peak_energies[0])
        colourplot(angles, peak_strengths, peak_strengths_error, (1 - optimised[0][0]) * full_dist2, t_angles, l_2, peak_energies[0])
        plt.plot(t_angles, sumdist, ls = 'dashed', color = 'xkcd:dark teal')

        #so you can't get the plots and simply paste them onto another set of axes, so we'll have to draw these axes again later. 
        dist_plotters = ['doublet', angles, peak_strengths, peak_strengths_error, optimised[0][0] * full_dist1, l_1, (1 - optimised[0][0]) * full_dist2, l_2, t_angles, peak_energies[0]]
        graphs.append(dist_plotters)
    elif l_select == 'n':
        pass
    else:
        colourplot(angles, peak_strengths, peak_strengths_error, n_dist_list[l_index][0], t_angles, l_select, peak_energies[0])
        #so you can't get the plots and simply paste them onto another set of axes, so we'll have to draw these axes again later. 
        dist_plotters = ['singlet',angles, peak_strengths, peak_strengths_error, n_dist_list[l_index][0], t_angles, l_select, peak_energies[0]]
        graphs.append(dist_plotters)
    
    plt.show()  

    if l_select == 'd':
        spectroscopicFactor1 = (n_dist_list[int(l_1)][0][0]/t_dist_list[int(l_1)][0]) * optimised[0][0]
        spectroscopicFactor2 = (n_dist_list[int(l_2)][0][0]/t_dist_list[int(l_2)][0]) * (1 - optimised[0][0])
        sfError1 = 0
        sfError2 = 0

        #add the first row
        row_dict = {'ENERGY': [peak_energy], 'L': [int(l_1)], 'SPECTROSCOPIC_FACTOR': [spectroscopicFactor1], 'ERROR': [sfError1]}
        row_df = pd.DataFrame(data = row_dict)   
        spectroscopic_df = spectroscopic_df.append(row_df)

        #and the second
        row_dict = {'ENERGY': [peak_energy], 'L': [int(l_2)], 'SPECTROSCOPIC_FACTOR': [spectroscopicFactor2], 'ERROR': [sfError2]}
        row_df = pd.DataFrame(data = row_dict)   
        spectroscopic_df = spectroscopic_df.append(row_df)

        print('The spectroscopic factor for the l = ', l_1, ' part of this state is:', spectroscopicFactor1)#, '\n\nThe theoretical distribution is: ', t_dist_list[l_index], '\n\nThe angles are: ', t_angles)
        print('The spectroscopic factor for the l = ', l_2, ' part of this state is:', spectroscopicFactor2)
    elif l_select == 'n':
        pass
    else:
        sfs = spectroscopic_finder(peak_strengths, peak_strengths_error, angles, t_dist_list[l_index], t_angles, int(l_select), n_dist_list[l_index][0][0]/t_dist_list[l_index][0], norm_errors_list[l_index])
        spectroscopicFactor = sfs[0]
        sfError = sfs[1]

        #add the row to the SF df that I'm making.
        row_dict = {'ENERGY': [peak_energy], 'L': [int(l_select)], 'SPECTROSCOPIC_FACTOR': [spectroscopicFactor], 'ERROR': [sfError]}
        row_df = pd.DataFrame(data = row_dict)   
        spectroscopic_df = spectroscopic_df.append(row_df)

        print('The spectroscopic factor for this state is:', spectroscopicFactor)#, '\n\nThe theoretical distribution is: ', t_dist_list[l_index], '\n\nThe angles are: ', t_angles)

    

    print('The peak cross-section for this state is:', max(peak_strengths))
    print('The position of this state is', peak_positions[0], '\n\n\n\n\n')


#print(graphs[0])
#print(peak_energies)

#get how many rows and columns there are
nplots = len(graphs)

npages = intdiv(nplots, 24)
os.chdir(current_directory)

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
            if graphs[i + pages * 24][0] == 'singlet':
                ax = colourplot(graphs[i + pages * 24][1], graphs[i + pages * 24][2],graphs[i + pages * 24][3],graphs[i + pages * 24][4], graphs[i + pages * 24][5], graphs[i + pages * 24][6],graphs[i + pages * 24][7] )
            elif graphs[i + pages * 24][0] == 'doublet':
                ax = colourplot_doublet(graphs[i + pages * 24][1], graphs[i + pages * 24][2],graphs[i + pages * 24][3],graphs[i + pages * 24][4], graphs[i + pages * 24][5], graphs[i + pages * 24][6],graphs[i + pages * 24][7],graphs[i + pages * 24][8],graphs[i + pages * 24][9])
            else:
                break                
        except:
            ax.axis('off')
            pass

        ax.set_xticks(np.arange(0,61,30))
        ax.set_xticks(np.arange(15,61,30), minor = True)
        #this code manually sets the y ticks. It is complicated and doesn't always work, and needs to be modified for doublets, but I'm keeping it around!
        '''    
        if max(graphs[i + pages * 24][2]) > max(graphs[i + pages * 24][4]):
            dist = graphs[i + pages * 24][1]
            errors = graphs[i + pages * 24][2]
            ax.set_yticks(ticklengthsetter(dist, errors))
            ax.set_yticks(ticklengthsetter(dist, errors, ticksno = 4), minor = True)

        else:
            dist = graphs[i + pages * 24][4]
            ax.set_yticks(ticklengthsetter(dist))
            ax.set_yticks(ticklengthsetter(dist, ticksno = 4), minor = True)
        '''
        ax.locator_params(tight=True, nbins=4)

        if i < 20 and i + pages * 24 < nplots -4:
            plt.setp(ax.get_xticklabels(), visible = False)

    fig.align_ylabels(axs=None)

    
    labelax.spines['top'].set_color('none')
    labelax.spines['bottom'].set_color('none')
    labelax.spines['left'].set_color('none')
    labelax.spines['right'].set_color('none')
    labelax.tick_params(labelcolor = 'w', top = False, bottom = False, left = False, right = False)
    labelax.set_xlabel('Angle (degrees)', size = 'xx-large', fontweight = 'bold')
    labelax.set_ylabel('Cross Section (mb/sr)', size = 'xx-large', fontweight = 'bold')

    fig.tight_layout()
    plt.show()

    #exit()

    if npages == 1:
        fig.savefig(reactionname + '_distribution.png')
    else:
        fig.savefig(reactionname + '_distribution_page_' + str(pages) + '.png')

    

spectroscopic_df = spectroscopic_df.reset_index(drop = True)
#print(spectroscopic_df)
os.chdir(analysis_code_directory + 'spectroscopic_factors/excel')
sfdf_filename = reactionname + '_spec_factors'
#print(sfdf_filename)

writer = pd.ExcelWriter(sfdf_filename + '.xlsx')
spectroscopic_df.to_excel(writer, 'Sheet1')
writer.save()
writer.close()

os.chdir(analysis_code_directory + 'spectroscopic_factors/pkl')
spectroscopic_df.to_pickle(sfdf_filename + '.pkl')
