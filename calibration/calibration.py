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

#I'm copying the database in functions, operating on that database, then returning the copy.
#Pandas doesn't like this as it thinks that I think that I'm working on the original.
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
#against position.
#For now, this doesn't include the uncertainties on the position or assigned energy.
#Therefore, it may not be completely accurate.
def calibrate(dat):
    
    #this finds all of the rows which have an energy assignment
    #it does this by checking whether the assigned energy on each row 'isfinite'
    #np.nan is not finite, and is what all of the energy row is defaulted to
    idx = np.isfinite(dat.EASSIGN)
    #this does the fit. It is the standard numpy polynmial fit, and uses least-squares regression
    #This doesn't take into account the uncertainties in energy, so it may need to be replaced
    #it will suffice for now.
    pl = np.polyfit(dat.POSITION[idx],dat.EASSIGN[idx],2)
    
    #makes points for line of best fit
    xp = np.linspace(min(dat.POSITION[idx]),max(dat.POSITION[idx]), max(dat.POSITION[idx]))
    
    #plots points and best fit line
    plt.plot(dat.POSITION[idx],dat.EASSIGN[idx],'o')
    plt.plot(xp,np.polyval(pl,xp),'g')

    #sets title and labels
    plt.title("Energy Calibration")
    plt.ylabel('Energy')
    plt.xlabel('Channel Number')

    plt.show()

    #now add an extra columb 'predicted energy'
    #this shows what the energy for each peak should be based on the current fit.
    predictedenergy = pl[0]*dat.POSITION**2 + pl[1]*dat.POSITION + pl[2]
    dat['PREDICTED_ENERGY'] = predictedenergy

    #finally, print off the important data after fitting so the user can evaluate what they 
    #want to do
    print(dat[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY']])

    return dat

def addrow(dat):

    #printing the data is tricky, because predicted energy is important.
    #If the predicted energy is present, print it. If not, print without it
    try:
        print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY']])
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
    savename = filename[0:-4] + '.pkl'
    
    #save it to a .pkl format
    savedata.to_pickle(savename)

    return finish

#This Function sets all of the predicted energies to all the final energies
def unassigned(dat):
    for i in range(dat.shape[0]):
        print(dat.EASSIGN[i])
        if np.isnan(dat.EASSIGN[i]) and dat.POSITION[i] != 0:
            dat.EASSIGN[i] = dat.PREDICTED_ENERGY[i]
            dat.sEASSIGN[i] = 0
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
    datadir = currentdir + '/peak_data/116Cd_d,p_peak_data/' #location of the data to load in

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
    data = pd.read_pickle(f)

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
            print(data[['POSITION','sPOSITION', 'AREA', 'sAREA','EASSIGN', 'sEASSIGN', 'PREDICTED_ENERGY']])
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



