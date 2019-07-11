#energy calibration.py
#the idea is to read in a set of peaks, and provide a nice interface to assign energies to those peaks
#as well as assiging energies, one should be able to draw a graph, and coefficients should be given.
#once this graph is drawn, one should be able to add more peaks.
#adding more peaks should be easier, sice a 'provisional' energy should be given for each peak to compare to ENSDF
#finally it should save the calibration coefficient and graph

#imports. os for changing directories and reading filenames, pandas for working with a database
#matplotlib.pyplot to plot calibration graphs and numpy for extra mathematical functionality
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.odr import *

#I'm copying the database in functions, operating on that database, then returning the copy.
#Pandas doesn't like this as it thinks that I think that I'm working on the original. It's a tool that is usually used for big data so copies aren't usually efficient in its intended use.
#Therefore it thinks I don't know what I'm actually doing, and complains at me.
#This setting tells it to shut up
pd.options.mode.chained_assignment = None


#A lot of the program is a loop once the data is loaded in, so if the user selects an option to
#quit, this flag will be switched to "true", which breaks the loop.
#S is a similar flag but dictates whether the user wants to save their calibration
finishflag = False


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FUNCTION DEFINITIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

#This function lists all of the files in the directory that is passed to it as an argument
#it also returns the number of files in said directory
def listfiles(directory):
    #set an index, since it's going to print out the index before each file so the user can
    #select which file to access more easily
    fileindex = 0
    #simple loop through the directory
    for roots, dirs, filenames in os.walk(directory):
        for f in filenames:
            fileindex = fileindex + 1
            #need to print the name of each file out so the user can select which one to open 
            print( str(fileindex) + ': ' + f)
    #returns the total number of files so the program can detect whether the
    #selected file number is valid (ie falls in between 0 and the number of files)
    return fileindex

#This function prompts the user to pick a file to open. It then checks the validity of the input.
#this technique is used to check the validity a few times in this module.
#the first step is to define the desired quantity, in this case the index of the file that
#the user wants to open. It initialises this quantity to None.
#A loop then begins, which asks the user to input the quantity.
#While the input is invalid, the loop continues.
#The input is checked using a try block. This try block tests the input by doing something
#that would throw an error if the argument was invalid. 
#If it passes this test, the value is returned by the function.
#If it throws an error, the 'except' block resets the quantity to None.
#The loop then continues.
def inputfileselector(fileno): #fileno is the number of files in the directory
    #Initialize the desired quantity
    fileselect = None
    while fileselect == None:
        #prompt the user to input their selection
        fileselect = input('Which file would you like to open?\n')
        try:
            
            #failure test
            #in this case, it checks to see if the input matches a file index
            
            #create a list of the possible file indices
            fileindices = list(range(1, fileno + 1))
            #filetest doesn't do annything, but if there's an error in the assignment, it
            #triggers the error
            filetest = fileindices.index(int(fileselect))

        except:
            fileselect = None
    return int(fileselect)

#a simple function which returns a filename.
#An index is passed to the function and the name of the file with that index is returned
def openfilefromlist(fileselect, directory):
    for roots, dirs, filenames in os.walk(directory):
        f = filenames[fileselect-1]
    return f

#Another selection function (see inputfileselector, line 50)
#This one prompts the user to select a row (in this case representing a peak)
#in a pandas dataframe.
def peakselector(DataFrame):
    #index of the desired peak
    peakselect = None
    while peakselect == None:
        #user selection prompt
        peakselect = input('Which peak would you like to assign an energy to?\n')
        try:
            #similar to inputfileselector, it checks whether the input matches a peak index 
            
            #DataFrame.shape[0] is the number of peaks
            peakindices = list(range(DataFrame.shape[0]))
            peaktest = peakindices.index(int(peakselect))

        except:
            #resets if there's an exception
            peakselect = None
    #returns the valid selection
    return int(peakselect)       

#Another selection function (see inputfileselector, line 50)
#This one asks for an energy assignment for a peak
def energyinput(peakindex, dat):
    #we want the energy and its uncertainty from this function
    energy = None
    error = None
    while energy == None and error == None:
        #user prompts
        energy = input('Input an energy in kev\n')
        error = input('Input an uncertainty on that energy, also in kev\n')
        try:
            #the test is to attempt to turn the input strings to floats.
            #Beware! Negative energies are allowed by the program, which is unphysical. 
            energy = float(energy)
            error = float(error)
        except ValueError:
            #gives an error message with this one
            print('At least one of the energy or uncertainty is not valid')
            energy = None
            error = None            
    #now we assign the valid energies and errors to the database, and return it
    dat.EASSIGN[peakindex] = energy
    dat.sEASSIGN[peakindex] = error
    return dat

#This function is the master function for deciding what the user wants to do.
#It is invoked once the user decides they don't want to assign more energies
#They will then get a list of options
#The function then calls the relevant function required to do the selected job
def jobselector(finishflag, dat):
    
    #This is the 'save' flag. If the user selects to save, this will be set to true
    #at the end of the loop, the finish flag will be set to 1 and the save function activated
    s = False

    #first step is to get the input
    decision = input('What would you like to do? Print table (\'p\'), set all energies to 0 (\'0\'), set remaining unassigned energies to predicted energies (\'u\'), do calibration and print graph (\'c\'), add a row(\'a\') quit (\'q\'), or save and finish (\'s\')\n')
    
    #now to deal with the inputs with if statemets
    
    #first we deal with the decision to print the table.
    if decision == 'p':
        print(dat)
        #simple enough! just prints the whole database

    #Next with them doing a calibration
    elif decision == 'c':
        #this, and the data, gets passed to a calibration function
        dat = calibrate(dat)
    
    #this is the quit option. It should just exit by setting the quit flag to 'True'    
    elif decision == 'q':
        finishflag = True

    #this is the 'save and quit' decision.
    #It sets a save flag instead of calling a function because it needs a lot of input arguments
    elif decision == 's':
        s = True

    #this sets all of the energy values to 0 in case one needs to do something quickly with uncalibrated data
    elif decision == '0':
        dat['EASSIGN'] = 0
        dat['sEASSIGN'] = 0
    
    #this sets all of the unassigned energy valuse to their predicted energies
    elif decision == 'u':
        try:
            dat = unassigned(dat)
        except KeyError: # if no calibration has been done, it 
            print('warning, no calibration done yet')

    elif decision == 'a':
        #this, and the data, gets passed to a function which adds an extra line to the file
        dat = addrow(dat)

    #everything else isn't on the options menu, and is therefore invalid.
    #In which case it should reset back to the start of the job selector
    #In this case this is done with recursion
    else:
        #error message
        print('Not a valid input')

        #now to reset to the start of this function. In this case it calls itself
        #when it does this it sets the return value of the inner function to the returning values
        #of this function, so when it does get out of the inner functions, the desired results
        #are returned by the outer function.
        recursivething = jobselector(finishflag, dat)
        dat = recursivething[0]
        finishflag = recursivething[1]
        s = recursivething[2]
    return [dat, finishflag, s]

#the actual calibration function. This does a quadratic fit on the assigned energy values
#against position. This takes into account x and y errors.
def calibrate(dat):
    
    #this finds all of the rows which have an energy assignment
    #it does this by checking whether the assigned energy on each row 'isfinite'
    #np.nan is not finite, and is what all of the energy row is defaulted to
    idx = np.isfinite(dat.EASSIGN)
    #do the fit using scipy's orthogonal distance regression
    #set our xdata and ydata first
    x = np.array(dat.POSITION[idx])
    y = np.array(dat.EASSIGN[idx])

    x_err = np.array(dat.sPOSITION[idx])
    y_err = np.array(dat.sEASSIGN[idx])

    y_err[y_err == 0] = 10**(-15)
            

    # Define a function (quadratic in our case) to fit the data with.
    def quad_func(p, x):
        return p[0]*x**2 + p[1] * x + p[2]
        #return p[0]*x**2 + p[1]

    # Create a model for fitting.
    quad_model = Model(quad_func)

    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, sx=x_err, sy=y_err)

    # Set up ODR with the model and data.
    odr = ODR(data, quad_model, beta0=[0.,1.,1.])

    # Run the regression.
    out = odr.run()

    # Use the in-built pprint method to give us results.
    out.pprint()
    '''Beta: [ 1.01781493  0.48498006]
    Beta Std Error: [ 0.00390799  0.03660941]
    Beta Covariance: [[ 0.00241322 -0.01420883]
    [-0.01420883  0.21177597]]
    Residual Variance: 0.00632861634898189
    Inverse Condition #: 0.4195196193536024
    Reason(s) for Halting:
    Sum of squares convergence'''

    x_fit = np.linspace(x[0], x[-1], 1000)
    y_fit = quad_func(out.beta, x_fit)

    plt.errorbar(x, y, xerr=x_err, yerr=y_err, linestyle='None', marker='x')
    plt.plot(x_fit, y_fit)


    #sets title and labels
    plt.title("Energy Calibration")
    plt.ylabel('Energy')
    plt.xlabel('Channel Number')

    plt.show()

    #now add an extra columb 'predicted energy'
    #this shows what the energy for each peak should be based on the current fit.
    predictedenergy = quad_func(out.beta, dat.POSITION)
    #error on this is: s_e^2 = (2ax + b)s_x^2 + x^4s_a^2 + x^2s_b^2 +s_c^2 + 2x^3cov_{ab} + 2x^2cov_{ac} + 2xcov_{bc}
    #apologies for the unreadable formula here but it is a long formula
    spredictedenergy = np.sqrt((2*out.beta[0]*dat.POSITION + out.beta[1])**2 * dat.sPOSITION**2 + dat.POSITION**4 * out.cov_beta[0][0] + dat.POSITION**2 * out.cov_beta[1][1] + out.cov_beta[2][2] + 2*dat.POSITION**3 * out.cov_beta[0][1] + 2*dat.POSITION**2 * out.cov_beta[0][2] + 2*dat.POSITION * out.cov_beta[1][2])
    dat['PREDICTED_ENERGY'] = predictedenergy
    dat['sPREDICTED_ENERGY'] = spredictedenergy

    #finally, print off the important data after fitting so the user can evaluate what they 
    #want to do
    print(dat[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY','sPREDICTED_ENERGY' ]])

    return dat

def addrow(dat):

    #printing the data is tricky, because predicted energy is important.
    #If the predicted energy is present, print it. If not, print without it
    try:
        print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY', 'sPREDICTED_ENERGY']])
    except KeyError:
        print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN']])
        

    peakix = float(input('Which peak do you want to add your new state before'))

    nPOSITION = 0
    nsPOSITION = 0
    nWIDTH = 0
    nsWIDTH = 0
    nHEIGHT = 0
    nsHEIGHT = 0
    nAREA = np.nan
    nsAREA = np.nan
    nCENTROID = 0
    nsCENTROID = 0
    nENERGY = float(input('What is the energy of this state'))
    nsENERGY = float(input('and what is its uncertainty'))
    nA = 0
    nsA = 0
    nB = 0
    nsB = 0
    nC = 0
    nsC = 0
    nR = 0
    nsR = 0
    nBETA = 0
    nsBETA = 0
    nS = 0
    nsS = 0



    line = pd.DataFrame({"POSITION": 0, "sPOSITION": 0, "AREA": 0, "sAREA": 0, "EASSIGN": nENERGY, "sEASSIGN": nsENERGY}, index=[peakix])
    dat2 = pd.concat([dat.ix[:peakix-1], line, dat.ix[peakix:]]).reset_index(drop=True)
    return dat2

#This is the saver function, which takes the database and saves the assigned energy peaks,
#along with other important data, to a file which is then accessed by the cross-section writer.
def saver(dat,finish, filename):
    #this is a save and quit option, so it will set the program to finish afterwards
    finish = True
    
    #now to save data
    
    #set the saving directory, then go there
    savedir = analysis_code_directory + 'spectrum_analysis/input'
    os.chdir(savedir)
   
    #need to chop down the database into the rows and columns that we want
    #first the rows, using the same technique as in calibrate()
    idx = np.isfinite(dat.EASSIGN)
    assigneddata = dat[idx].copy()
    
    #next we chop out the predicted energy column
    savedata = assigneddata[['POSITION', 'sPOSITION', 'WIDTH', 'sWIDTH', 'HEIGHT' ,'sHEIGHT', 'AREA', 'sAREA', 'CENTROID', 'sCENTROID', 'ENERGY', 'sENERGY', 'A', 'sA', 'B', 'sB', 'C', 'sC', 'R', 'sR', 'BETA', 'sBETA', 'S', 'sS', 'EASSIGN', 'sEASSIGN' ]].copy()

    savedata = savedata.reset_index(drop = True) 
    #set the filename by chopping off the .txt and puting .pkl instead
    savename = filename[0:-4] + '.txt'
    

    #save it to a .pkl format
    #savedata.to_pickle(savename)
    savedata.to_csv(savename, sep = ' ')

    return finish

#This Function sets all of the predicted energies to all the final energies
def unassigned(dat):
    for i in range(dat.shape[0]):
        print(dat.EASSIGN[i])
        #want to assign all of the unassigned energies with what the fit says they should be
        if np.isnan(dat.EASSIGN[i]) and dat.POSITION[i] != 0:
            dat.EASSIGN[i] = dat.PREDICTED_ENERGY[i]
            dat.sEASSIGN[i] = dat.sPREDICTED_ENERGY[i]
        #if the fit is more accurate than previous results, replace the previous result with the new one.
        if dat.sEASSIGN[i] > dat.sPREDICTED_ENERGY[i]:
            dat.EASSIGN[i] = dat.PREDICTED_ENERGY[i]
            dat.sEASSIGN[i] = dat.sPREDICTED_ENERGY[i]
    return dat

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN PROGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

#firstly need to tell the program where it is and where to load things from
currentdir = os.getcwd()
analysis_code_directory = currentdir[0:-11] #this allows the program to work via relative paths

#need to set up a choice of if they want to do new calibration or an old one
newOrLoad = input('Would you like to start a new calibration, or load an existing one? n/l')

if newOrLoad == 'n':
    loadmethod = "new"
    datadir = currentdir + '/peak_data/116Cd_h,a_peak_data/' #location of the data to load in

    #first get the number of files in the directory
    numberoffiles = listfiles(datadir)

    #then get the user to select which file to open, and open it
    fileselect = inputfileselector(numberoffiles)
    f = open(datadir + openfilefromlist(fileselect, datadir))   

    #read in the data from that file
    data = pd.read_table(f, sep = ' ')
    f.close()
    #add columns for assigned energy and its error
    data['EASSIGN'] = np.nan
    data['sEASSIGN'] = np.nan

elif newOrLoad == 'l':
    loadmethod = "load"
    loaddir = analysis_code_directory + 'spectrum_analysis/input/'
    os.chdir(loaddir)

    #first get the number of files in the directory
    numberoffiles = listfiles(loaddir)

    #then get the user to select which file to open, and open it
    fileselect = inputfileselector(numberoffiles)
    f = loaddir + openfilefromlist(fileselect, loaddir)
    
    #read in data, takes filepath as its argument
    try:
        data = pd.read_pickle(f)
    except:
        data = pd.read_table(f, sep = ' ')

#It's early in the program at this stage so I'll just quit if there's invalid input
else:
    print('Error, no valid option selected, quitting program.')
    exit()

#set a loop. this is because the user will want to add multiple assignments
#do multiple calibrations etc.
while finishflag == False:
    #ask if they'd like to assign a state
    assignflag = input('Would you like to assign an energy to a state?y/n\n')
    
    #if so, show them the states with important data (position, predicted energy, strength,
    #assigned energy), then ask them to make the assignment
    if assignflag == 'y':
        
        #printing the data is tricky, because predicted energy is important.
        #If the predicted energy is present, print it. If not, print without it
        try:
            print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY', 'sPREDICTED_ENERGY']])
        except KeyError:
            print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN']])
        
        #now get the user to pick a peak
        peakindex = peakselector(data)
        #and assign an energy to it
        data = energyinput(peakindex, data)

    #if they don't want to assign energy, give the user a choice of what they want to do
    #print the table, calibrate, save and quit, or quit.
    elif assignflag == 'n':
        #run the function that does this. It returns a copy of the updated data
        #it also returns finish and save flags
        jobselected = jobselector(finishflag, data)
        data = jobselected[0]
        finishflag = jobselected[1]

        #if th =e save flag is true, run the save program.
        #this sets a True finish flag as well
        if jobselected[2] == True:
            #save over original name if loaded, save newly if a new calibration
            if loadmethod == "new":
                finishflag = saver(data,finishflag,openfilefromlist(fileselect, datadir)) 
            else:
                finishflag = saver(data,finishflag,f)
    else:
        #if they don't put yes or no, it prints an error message and continues the loop 
        print('Not a valid input')



