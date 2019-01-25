#extraction.py
#gets the raw data and converts it to a form that is easier to read by people and other programs, then saves it in another directory

import os


#functions
def betweenbrackets(astring):
    index1 = astring.index('(')
    index2 = astring.index(')')
    return astring[index1+1:index2]   

def beforebrackets(astring):
    bracketindex = astring.index('(')
    return astring[0:bracketindex]

def errorfinder(astring):

    index1 = astring.index('(')
    index2 = astring.index(')')

    try:

        decimalindex = astring.index('.')
        decplac = len(astring[decimalindex+1:index1])
        afterdp  = float(betweenbrackets(astring))
       
        error = str(afterdp / (10 **decplac))
    
    except ValueError:

        error = astring[index1+1:index2]

    return error
    
#set the current directory and analysis code directories
currentdir = os.getcwd()
analysisdir = currentdir[:-15]   

#set target and reaction
target = "116Cd"
reaction = "d,p"

#fset the directory to convert and publish to
indir = "%speak_extraction/raw_data_%s_%s_new"%(analysisdir,target,reaction)
odir = "%scalibration/peak_data"%(analysisdir)

#now go to output directory, and make a folder there to put the sorted data to
os.chdir(odir)
os.mkdir("%s_%s_peak_data_new"%(target,reaction))
odir = odir + "/%s_%s_peak_data_new"%(target,reaction)

os.chdir(indir)

#now we want it to loop through the files in the directory
for root, dirs, filenames in os.walk(indir):
    for file in filenames:
        print(file)
        #get all the text in a file on a string so it can be sifted through
        f = open(file)
        filestring = f.read()
        f.close()
        filestring = filestring.split()

        #now go to output directory, and make a folder there to put the sorted data to
        os.chdir(odir)
        ofilename = file[0:-8] + ".txt"
        with open (ofilename, 'w') as ofile:
            
            header_line = 'POSITION sPOSITION WIDTH sWIDTH HEIGHT sHEIGHT AREA sAREA CENTROID sCENTROID ENERGY sENERGY A sA B sB C sC R sR BETA sBETA S sS\n'
            ofile.write(header_line)
           
            #loop over all of the data files
            for i in range(len(filestring)):
                #we want to isolate each individual fit to extract its properties.
                #The delimiter for fits is the word 'Background:'
                #because it is the first word in the gf3 output
                #so we get the index i equaling the start of a fit and the index j being the end
                if filestring[i] == "Background:":
                    j = i+1
                    try:
                        #use a loop to find the index j which is the next instance of 'Background'
                        while filestring[j] != "Background:":
                            j = j + 1
                    except IndexError: 
                        #of course if j = i+1 and goes the length of the file, j will sometimes be out of range
                        pass
                
                    #print(filestring[i:j])       
                    #assign the snippet corresponding to the fit to a string
                    fit = filestring[i:j]
                
                    #there are 3 cases, firstly where there is no peak, Background: -1
                    #secondly where there is a sum and position but no other information
                    #thirdly where there are peak fits
                    try: #this is to catch the case where at the end of the file there's just Background 
                        if fit[1] == '-1':
                            pass
                    except IndexError:
                        continue
                    
                    #in this case, there is no peak detected where there is a peak here at other angles
                    #in this case we need to write a null line in the output file so the other
                    #analysis knows there is no peak for this angle 
                    if fit[1] == '-1': 
                        ofile.write("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")
                    elif fit[1] == 'POS':
                        POSITION = beforebrackets(fit[2])
                        sPOSITION = betweenbrackets(fit[2])
                        AREA = beforebrackets(fit[4])
                        sAREA = betweenbrackets(fit[4])
                        ofile.write("%s %s 0 0 0 0 %s %s 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"%(POSITION, sPOSITION, AREA, sAREA))

                    elif fit[1] == 'A':
                        #this is the case where there is a fit
                        #there may be multiple peaks for a fit
                        #first extract the general paramerters for a fit
                        #print(fit)
                        a = beforebrackets(fit[3])
                        sa = errorfinder(fit[3])
                        b = beforebrackets(fit[6])
                        sb = errorfinder(fit[6])
                        c = beforebrackets(fit[9])
                        sc = errorfinder(fit[9])
                        r = beforebrackets(fit[13])
                        sr = errorfinder(fit[13])
                        beta = beforebrackets(fit[16])
                        sbeta = errorfinder(fit[16])
                        s = beforebrackets(fit[16])
                        ss = errorfinder(fit[16])
                            
                        #now need to loop over the individual peaks in the fit
                        #every numerical quantity should have an error associated with it
                        #except the peak number
                        peakno = 1
                        for elements in range(len(fit)):
                            if fit[elements] == str(peakno):
                                #print(elements)
                                position = beforebrackets(fit[elements + 1])
                                sposition = errorfinder(fit[elements + 1])
                                width = beforebrackets(fit[elements + 2])
                                swidth = errorfinder(fit[elements + 2])
                                height = beforebrackets(fit[elements + 3])
                                sheight = errorfinder(fit[elements + 3]) 
                                area = beforebrackets(fit[elements + 4])
                                sarea = errorfinder(fit[elements + 4])
                                centroid = beforebrackets(fit[elements + 5])
                                scentroid = errorfinder(fit[elements + 5])
                                energy = beforebrackets(fit[elements + 6])
                                senergy = errorfinder(fit[elements + 6])

                                ofile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n"%(position, sposition, width, swidth, height, sheight, area, sarea, centroid, scentroid, energy, senergy, a, sa, b, sb, c, sc, r, sr, beta, sbeta, s, ss))


                                peakno = peakno + 1
                        
                        
                    
            os.chdir(indir)

