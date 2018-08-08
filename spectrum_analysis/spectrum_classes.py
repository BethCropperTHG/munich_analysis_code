import math


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~DEFINE SPECTRUM AND PEAK CLASSES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#Define a class which contains infomation about a peak
#this needs to come first since spectrum class should have functions that work on peaks
class peak:
	#constructor
	def __init__(self,POSITION,P_ERR,AREA,A_ERR,ENERGY,E_ERR):

		self.position = POSITION
		self.p_err = P_ERR
		self.area = AREA
		self.a_err = A_ERR
		self.energy = ENERGY
		self.e_err = E_ERR

#define a spectrum class which has all of the spectrum data on
class spectrum:
	#constructor
	def __init__(self,name,TARGET_THICKNESS,THICKNESS_UNCERTAINTY,OMEGA,THETA,S1,S3,CURRENT_OFFSET,CURRENT_GRADIENT,RUN_TIME,FS,ZBEAM, ZTARGET,ATARGET,CH1,e_Ch1,Rest,e_R,ISOTOPIC_PURITY):

		#stick the constructor paramaters onto the object variables
		self.name = name
		self.target_thickness = TARGET_THICKNESS
		self.thickness_uncertainty = THICKNESS_UNCERTAINTY
		self.omega = OMEGA
		self.theta = THETA
		self.s1 = S1
		self.s3 = S3
		self.current_offset = CURRENT_OFFSET #
		self.current_gradient = CURRENT_GRADIENT
		self.run_time = RUN_TIME
		self.fs = FS
		self.zbeam = ZBEAM
		self.ztarget = ZTARGET
		self.atarget = ATARGET
		self.ch1=CH1
		self.e_ch1 = e_Ch1
		self.rest = Rest
		self.e_r = e_R
		self.isotopic_purity = ISOTOPIC_PURITY

		#derived variables
		self.beam_offset = self.current_offset * self.run_time
		self.bic_corrected = ((((self.s1 - self.s3)*self.fs*10**9)/(1000*self.run_time))+self.current_offset*1.6e-19*10**9)*self.run_time*0.5
		#self.no_beam = ((self.s1-self.s3)*self.fs)/(self.zbeam*1.61e-19*1000)
		self.no_beam = (self.bic_corrected*self.fs)/(self.zbeam*1.61e-19*1000)
		
		self.sno_beam = ((math.sqrt(self.s1)+math.sqrt(self.s3))*self.fs)/(self.zbeam*1.61e-19*1000)
		self.no_target = self.target_thickness/(self.atarget*1.661e-18) #per millibarn
		self.snotarget = self.thickness_uncertainty/(self.atarget*1.661e-18)
		self.efficiency = self.rest/(self.ch1+self.rest)
		self.sefficiency = self.efficiency*math.sqrt((self.e_r/self.rest)**2+((self.e_ch1**2+self.e_r**2)/(self.ch1+self.rest))**2)
		self.xs_sys_error = math.sqrt((self.sefficiency/self.efficiency)**2+(self.snotarget/self.no_target)**2+(self.sno_beam/self.no_beam)**2)

		#declared but unused variables
		self.xsection = None
		self.sxsection = None

		#needs a list of peaks assigned to the class
		self.peaks = []

		'''
		print(self.name)
		print(self.s1-self.s3)
		print(self.bic_corrected)
		print(self.fs)
		print(self.efficiency)
		print(self.no_beam)
		print('\n')
		print(self.run_time)
		print(self.current_offset)
		print('\n')
		'''
	#needs a function to operate on peaks in the spectrum to find cross section


	def peak_xsection(self,i):
		self.peaks[i].xsection = (1e27*self.peaks[i].area)/(self.no_target*self.no_beam*(self.omega/1000)*self.isotopic_purity*self.efficiency)

		self.peaks[i].sxsection = (1e27*self.peaks[i].a_err)/(self.no_target*self.no_beam*(self.omega/1000)*self.isotopic_purity*self.efficiency)   
		'''
		if self.theta < 11:
			self.peaks[i].xsection = self.peaks[i].xsection * 2
			self.peaks[i].sxsection = self.peaks[i].sxsection * 2
		
		if i == 0:
			print(self.name)
			print(self.peaks[i].energy)
			print(self.peaks[i].area)
			print(self.no_target)
			print(self.no_beam)
			print(self.omega)
			print(self.peaks[i].xsection)
			print('\n')
		else:
			pass
		'''
		return;



















