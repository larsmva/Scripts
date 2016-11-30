#!/usr/bin/env python
# -*- coding: utf-8 -*-
from dolfin import *

"""
DICOM header : http://dicomlookup.com/lookup.asp?sw=Ttable&q=C.8-4


Comments: 

	There exsist different cut-off concentration for different tissue/fluid.

	The cut-off concentration is defined as 
		
			(math.log( (TR*r1+TE*r2)/(TE*r2) ) -TR*R01)/(TR*r1) --> Check

	and is the value for the signal derivative of the conctration is zero.
	Thus the maximum signal intensity possible is obtained at the cut-off 
	concentration, at larger concentration the signal intensities decays.
	

To do :

	define affine transform -> same for all 
"""


def mric_solver_SE(k1,k2,G,Is ,max_num=3000 , tol=1.0e-9, cut_off=True): # TO DO cut_off of I
	"""
	Description:
		    Solver for F(c) =  ..... 

		    The derivative of F(c) is modified 
	Comments:
		There are different formulas for different scan sequences. This mehtod is 
		for spin echo (SE). Type of scan sequence can be read from dicom header.
	==============	 ==============================================================
	Argument                Explanation
	==============	 ==============================================================
	k1		 Defined as TR*r1
	k2		 defind as TE*r2
	G		 defined exp( -TR*R1)
	I		 The relative intensities 

	"""
	import numpy as np
	I=np.copy(Is) # avoid changing the array	
	c_max = -np.log( k2/(G*k1+G*k2)) /k1 
	c_ = c_max*((I-1.)/I )  # initial guess


	def F(c):
		return (1.-G*np.exp(-k1*c))*np.exp(-k2*c) -I*(1.-G) 
	def H(c):
		return - F(c)/(k2-(k1+k2)*G*np.exp(-k1*c)) # add /np.exp(-k2*c) ?? to get ->F(c)/F'(c)
	def RI(c):
		return (1.-G*np.exp(-k1*c))*np.exp(-k2*c)/(1.-G) 

	# CUT_OFF
	if cut_off==True:
		max_I = RI(c_max) 
		I[I>max_I] = max_I

	num_iter=0
	e = 1.0
	while e > tol and num_iter < max_num :
		c = c_ - H(c_)
		num_iter+=1
		e = abs(c_ - c).max()  # ensures that all values must converge Maybe RI(c) -I  is better ? 
		c_=c
		if num_iter > max_num:
			raise ValueError('Convergence error')

	return c




def mric_brentq(k1,k2,G,theta,I,max_num=3000 , tol=1.0e-7): #TO DO improve, implement with c++ boost
	
	def Ir_GRE(c) : 
		return (1.-aux0*G)*(1.-G*np.exp(-k1*c))*np.exp(-k2*c)/((1.-aux0*G*np.exp(-k1*c))*(1.-G))

	import numpy as np
	
	import scipy
	aux0 = np.cos(np.deg2rad(theta)) 
	aux1= 1.-aux0
	aux2 = 1.+aux0
	x = ( k1*aux1 + k2*aux2  + np.sqrt(k1**2*aux1**2 + 2*k1*k2*aux1*aux2 + k2**2*aux1**2))/(2*k2) # if /aux0
	c_max = np.log(x/G)/k1  	



	
	Im = Ir_GRE(c_max) 

	if I>Im:
	   I = Im
	
	if I<1.0:
	   I=1.0

	def F(c) : 
		return (1.-G*np.exp(-k1*c))*np.exp(-k2*c)/(1.-aux0*G*np.exp(-k1*c))  -I*(1.-G)/(1.-aux0*G) 
	

	
	c =  scipy.optimize.brentq(F,0.0,c_max)
	
	
	return c
	


def mric_GRE_solver(k1,k2,G,theta,I,max_num=3000 , tol=1.0e-7):
	return False

	

	
	
	




def get_header_info(path2dicom):
	"""
	Reads the listed info from the dicom-file:
		TE  = echo time
		TR  = Repetition time
		magnetric strength
		scanning sequence  
	""" 
	import dicom # REQUIRE DICOM
	ds = dicom.read_file(path2dicom)
	TE = float(ds.EchoTime)/1000 		#ms->s
	TR = float(ds.RepetitionTime)/1000      #ms->s
	mfs = float(ds.MagneticFieldStrength)
	ss = str(ds.ScanningSequence)	
	if ss=="GR":	
		fa = float(ds.FlipAngle)
	else : 
		fa = 0
	
	r1,r2,R01,R02 = values_from_database(mfs)
	
	return TR,TE,r1,r2,R01,R02,ss,fa

def values_from_database(magnetic_strength,tissue="white"): 
	# r1 and r2 taken from MICAD : Molecular Imaging and Contrast Agent Database  DIMENSION ARE l/(mmol s)
	# R01 and R02 taken from : T1, T2 relaxation and magnetization transfer in tissue at 3T , Stanisz, Greg J. et al. * 1/s
	
	v = {"white":{1.5:[4.7,6.8,1.13,13.9],3.0:[ 3.6 , 6.3,0.92, 14.5]},"grey":{1.5:[4.7,6.8,0.89,10.5],3.0:[3.6,6.3,0.55,10.1] } }

	#r1     r2       R01  	R02	
	return v[tissue][magnetic_strength]






def neighboor_values( data,i,j,k):
	import numpy as np
	
	temp = data[i-2:i+3,j-2:j+3,k-2:k+3].reshape(1,-1)[0]
	v = np.sort(temp)
	return sum(v[-6:-1])/len(v[-6:-1])

		

def get_MRI_values(path2mri_file,function_space):
		"""
		==============	 ==============================================================
		Argument                Explanation
		==============	 ==============================================================
		path2_mri_file	The path to the MRI file with extension mgz. typically orig.mgz
		function_space	The function space of the mesh.		
		mesh 	 	The mesh of the brain	 
		function	Instead of a return value, update the function inside. 	
		offset 		The translation of the origin of the mesh , need to be np.array
		"""
		import nibabel
		from nibabel.affines import apply_affine
		import numpy.linalg as npl
		import numpy as np

		img= nibabel.load(path2mri_file) 
		inv_aff = npl.inv ( img.get_header().get_vox2ras_tkr() )
		data = img.get_data()
	

		xyz = function_space.tabulate_dof_coordinates().reshape((function_space.dim(),-1))

		i,j,k = apply_affine(inv_aff,xyz).T
	
		i= map(round,i) 
		j= map(round,j) 
		k= map(round,k) 
		
		values = np.zeros(len(i))
		for n in range(len(i)):
			values[n]= neighboor_values(data,i[n],j[n],k[n]) 



		return values







