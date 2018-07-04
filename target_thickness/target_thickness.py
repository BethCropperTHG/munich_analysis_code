import numpy as np
import pandas as pd
import math
import os

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET GLOBAL DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#define a directory in which data is
data_directory = os.getcwd() + '/'

#open the experiment properties file, and extract all of the general properties of the system
with open('%sEXPT_PROPERTIES.txt'%data_directory) as prop:
    propertiesstring = prop.readlines()[1]
    properties = propertiesstring.split()

#to make this more easy to read I'll store these in individual variables
theta_degrees = float(properties[0])
fc_scale = float(properties[1])
solid_angle_msr = float(properties[2])
e_mev = float(properties[3])
z_beam = float(properties[4])

def degreestoradians(angle): return (angle * 2 * 3.141592654)/360

#some derived definitions
theta = degreestoradians(theta_degrees)
solid_angle = solid_angle_msr/1000

"""
~~~~~~~~~~~~~~~~~~~MAKE TARGET CLASS TO CALCULATE TARGET THICKNESS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

class target:
    #initialise with necessary properties: A,Z,S1,S3,Yield, and Yield uncertainty
    #annoyingly, yield is a thing in python so we'll capitalise our Yield
    def __init__(self,A,Z,element,S1,S3,Yield,yield_uncertainty, isotopic_purity, efficiency):

        #first, define obvious properties
        self.A = A
        self.Z = Z
        self.S1 = S1
        self.element = element
        self.S3 = S3
        self.Yield = Yield
        self.yield_uncertainty = yield_uncertainty
        self.isotopic_purity = isotopic_purity
        self.efficiency = efficiency

        #now need calculated ones
        #rutherford cross section, nothing too bad here, using nice formula with alpha and h_bar*c in it
        self.rutherford = 1e-30*((self.Z*z_beam*(1/137)*197.327)**2)/(4*e_mev*np.sin(theta/2)**2)**2

        #with #beam, factor of 1000 to fc scale, weird unit changes too
        self.no_beam = ((self.S1-self.S3)*fc_scale)/(z_beam*1.61e-19*1000)

        #error in number in the beam
        self.sigma_no_beam = ((math.sqrt(self.S1)+math.sqrt(self.S3))*fc_scale)/(z_beam*1.61e-19*1000)

        #this is the cross section formula, dsigma/domega = y/#target#beam*solid angle*efficiency
        self.no_target = self.isotopic_purity * self.Yield/(self.rutherford*self.no_beam*solid_angle*efficiency)

        #this is basically unit conversion from #nuclei/m**2 to micrograms/cm**3
        self.target_thickness = (self.no_target/10000)*self.A*1.661e-27*1000*1000000

        #time for error, basically fractional uncertainties with the yield uncertainty.
        #this is added to the beam number uncertainty (scaler uncertainties)
        #these have counting statistics, sigma(s1) = sqrt(s1)
        self.thickness_uncertainty = math.sqrt((self.yield_uncertainty/self.Yield)**2 + (self.sigma_no_beam/self.no_beam)**2)*self.target_thickness

"""
~~~~~~~~~~~~~~~READ IN DATA, CALCULATE TARGET THICKNESSES, THEN SAVE EVERYTHING~~~~~~~~~~~~~~~
"""

#open file and stick the header line there with the columns I want
f = open('TARGET_PROPERTIES.txt', 'w')
f.write('TARGET PURITY TARGET_THICKNESS UNCERTAINTY\n')
f.close()

#create a pandas table based on the RUN_PROPERTIES file, i.e the results for each target
#pandas is nice and automatically titles the columns with the header
targets_df = pd.read_table('%sRUN_PROPERTIES.txt'%data_directory, sep = ' ')


#this 'shape' is a list of the length and width of the table
#looping doesn't like accessing class data so we save it before the loop
shape = targets_df.shape

#loop over the rows
for i in range (0,shape[0]):
    
    #read in data from the row
    t_A = targets_df.A[i]
    t_Z = targets_df.Z[i]
    t_element = targets_df.ELEMENT[i]
    t_S1 = targets_df.S1[i]
    t_S3 = targets_df.S3[i]
    t_YIELD = targets_df.YIELD[i]
    t_YIELD_UNCERTAINTY = targets_df.YIELD_UNCERTAINTY[i]
    t_ISOTOPIC_PURITY = targets_df.ISOTOPIC_PURITY[i]
    t_EFFICIENCY = targets_df.EFFICIENCY[i]

    #create a 'target' object which will automatically calculate target thickness and uncertainties
    t_object = target(t_A,t_Z,t_element,t_S1,t_S3,t_YIELD,t_YIELD_UNCERTAINTY,t_ISOTOPIC_PURITY,t_EFFICIENCY)
    
    #write the relevant data from that object to a file (element, isotopic purity of target, target thickness, uncertainty)
    f = open('TARGET_PROPERTIES.txt', 'a')
    f.write('%s%s %s %s %s\n'%(t_A,t_element,t_ISOTOPIC_PURITY,t_object.target_thickness,t_object.thickness_uncertainty))

    f.close()























