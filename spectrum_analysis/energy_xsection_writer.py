import numpy as np
import pandas as pd
import os
import spectrum_classes as sp


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET SPECTRUM DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#define the directory where we are
current_directory = os.getcwd() + '/'

#get empty list to stick spectrum objects on
spectra = []


os.chdir('run_properties')

#loop over the files in the directory and read spectrum data to an individual spectrum object
#then we want to stick those objects on a list to be accessed later
indir = '%srun_properties'%(current_directory)
#I don't know what this does, I robbed it off stack overflow
#seems to loop over files in a directory
for root, dirs, filenames in os.walk(indir):
	#apparently f is a string of a filename, and filenames is a list/tuple of filenames
	for f in filenames:
		#print(filenames)
		#open the experiment properties file, and extract all of the general properties of the system
		with open(f) as prop:
			
			#split the file into lines, then split that into numbers
			lines = list(prop)
			print('The lines of filename ', f, ' are:')			
			print(lines, '\n')
			
			propertiesstring = lines[1]
			#print(propertiesstring)
			properties = propertiesstring.split(',')   
			print('The list of properties is then therefore:', properties, '\n')
			#print(properties[3][1:-1])
			#to make this more easy to read I'll store these in individual variables
			try:
				TARGET_THICKNESS = float(properties[0])	
				THICKNESS_UNCERTAINTY = float(properties[1])
				OMEGA = float(properties[2])
				THETA = float(properties[3])
				S1 = float(properties[4])
				S3 = float(properties[5])
				CURRENT_OFFSET = float(properties[6])
				CURRENT_GRADIENT = float(properties[7])
				RUN_TIME = float(properties[8])
				FS = float(properties[9])
				ZBEAM = float(properties[10])
				ZTARGET = float(properties[11])
				ATARGET = float(properties[12])
				CH1 = float(properties[13])
				e_CH1 = float(properties[14])
				REST = float(properties[15])
				e_R = float(properties[16])
				ISOTOPIC_PURITY = float(properties[17][0:-1])#get rid of newline character

			except ValueError:

				TARGET_THICKNESS = float(properties[0][1:-1]) #need these [1:-1] because the way they split includes some quotation marks	
				THICKNESS_UNCERTAINTY = float(properties[1][1:-1])
				OMEGA = float(properties[2][1:-1])
				THETA = float(properties[3][1:-1])
				S1 = float(properties[4][1:-1])
				S3 = float(properties[5][1:-1])
				CURRENT_OFFSET = float(properties[6][1:-1])
				CURRENT_GRADIENT = float(properties[7][1:-1])
				RUN_TIME = float(properties[8][1:-1])
				FS = float(properties[9][1:-1])
				ZBEAM = float(properties[10][1:-1])
				ZTARGET = float(properties[11][1:-1])
				ATARGET = float(properties[12][1:-1])
				CH1 = float(properties[13][1:-1])
				e_CH1 = float(properties[14][1:-1])
				REST = float(properties[15][1:-1])
				e_R = float(properties[16][1:-1])
				ISOTOPIC_PURITY = float(properties[17][1:-2])#get rid of newline character


			#create spectrum object with all the variables I just made dumped into the constructor
			spectrum_name = f[0:-14]
			this_spectrum = sp.spectrum(spectrum_name,TARGET_THICKNESS,THICKNESS_UNCERTAINTY,OMEGA,THETA,S1,S3,CURRENT_OFFSET,CURRENT_GRADIENT,RUN_TIME,FS,ZBEAM, ZTARGET, ATARGET, CH1,e_CH1,REST,e_R,ISOTOPIC_PURITY)

			#append spectrum to list
			spectra.append(this_spectrum)


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~READ IN PEAK DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
input_dir = '%s/input'%current_directory
os.chdir(input_dir)


#yay, nested loops! loops over the folder, input file names, and spectra to match a input file up to a 
#spectrum. then we create a pandas dataframe for the peak data and read that onto a peak
for root, dirs, filenames in os.walk(input_dir):
    for f in filenames:
        for i in range(len(spectra)):
            print(spectra[i].name + '    ' + f[0:29])
            if spectra[i].name == f[0:-9]:
                print('The name of the spectrum is: ', spectra[i].name)
				#create a pandas table based on the file
                try:
                    peaks_df = pd.read_table(f, sep = ' ')
                    print('The file is a .txt file')
                except:
                    peaks_df = pd.read_pickle(f)
                    print('The file is a .pkl file')
                #print('\n\n\n', peaks_df, '\n\n\n')
                #this 'shape' is a list of the length and width of the table
                #looping doesn't like accessing class data so we save it before the loop
                shape = peaks_df.shape

                #loop over the rows
                for j in range (shape[0]):
					
                    #read in data from the row
                    POSITION = peaks_df.POSITION[j]
                    sPOSITION = peaks_df.sPOSITION[j]
                    AREA = peaks_df.AREA[j]
                    sAREA = peaks_df.sAREA[j]
                    ENERGY = peaks_df.EASSIGN[j]
                    sENERGY = peaks_df.sEASSIGN[j]


                                        

                    #create a peak object to stick on the peaks list in the spectrum
                    peak_temp = sp.peak(POSITION,sPOSITION,AREA,sAREA,ENERGY,sENERGY)


                    #stick the saved peak on the peaks list in the current spectrum
                    spectra[i].peaks.append(peak_temp)
                    
                    #work out energies and cross sections for the peak
					
                    spectra[i].peak_xsection(j)
                    #print('The cross section for the ', j, 'th peak of the ', i, 'th spectrum is: ', spectra[i].peaks[j].xsection, '\n')

#need to save the spectrum to a file, 1 file per spectrum, 1 row per peak

#change to output directory
output_dir = '%s/output'%current_directory
os.chdir(output_dir)

for i in range(len(spectra)):

	#create a file for a spectrum
	with open ('%speaks.txt'%spectra[i].name, 'w') as f:

		#write a header line to the file
		header_line = 'PEAK_POSITION PEAK_ENERGY CROSS_SECTION ERROR\n'
		f.write(header_line)
		print('Writing header file to file to ', f)

		#now want to loop over the peaks and add peak information
		for j in range(len(spectra[i].peaks)):
			print('The position of this peak is at:',spectra[i].peaks[j].position)
			f.write('%s %s %s %s\n'%(spectra[i].peaks[j].position, spectra[i].peaks[j].energy, spectra[i].peaks[j].xsection, spectra[i].peaks[j].sxsection,))
