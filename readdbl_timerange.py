#Usage : Call the function as readbinary(directory_to_dblfiles,Resolution,Datatype_flag,tbeginnging,tend)
# Flags : 1-Density , 2-vx1 , 3-vx2 , 4-vx3 , 5-Tgas,6-ionx , 7-iony , 8-prs , 9-ueuv , 10-urec , 11 -fun,12-bx1,13-bx2,14-bx3, 15-phi


from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import ReferenceUnits as ref

def readsinglefile(filedirectory,timestep,N,field,code_units=False) : 
	"""
	File read should have path filedirectory/field.timestep.dbl
		filedirectory: string 
			directory in which files reside
		timestep : integer 
			timestep of file
		N : integer
			Resolution
		field: string
			Data field to read in
		code_units : Boolean 
			Flag to return in code units
	"""

	if(field == 'rho') : field = "rho."
	elif(field == 'vx1') : field = "vx1."
	elif(field == 'vx2') : field = "vx2."
	elif(field == 'vx3') : field = "vx3."
	elif(field == 'Tgas') : field = "Tgas."
	elif(field == 'ionx') : field = "ionx."
	elif(field == 'iony') : field = "iony."
	elif(field == 'prs') : field = "prs."
	elif(field == 'ueuv') : field ="ueuv."
	elif(field == 'urec') : field = "urec."
	elif(field == 'fun') : field = "fun."
	elif(field == 'bx1') : field = "bx1."
	elif(field == 'bx2') : field = "bx2."
	elif(field == 'bx3') : field = "bx3."
	elif(field == 'phi') : field = "phi."
	else : raise ValueError('The field flag does not exist')

	filedirectory = os.path.abspath(filedirectory) + '/'

	filename = filedirectory+field+"%04d"%timestep+".dbl"

	try:
			file = open(filename,'rb')
	except:
			print ('File %s Not found'%(filename))
			raise SystemExit

	#
	# Read binary file
	#

	file.seek(0,2)	
	eof = file.tell()
	file.seek(0,0)

	shape = (N,N,N)
	count = prod(shape)

	data = fromfile(file,dtype=float,count=count)

	if(code_units):
		unit_factor = 1.0
	else:
		if(field == 'rho.') : unit_factor = ref.unit_Density
		elif(field == 'vx1.') : unit_factor = ref.unit_Velocity
		elif(field == 'vx2.') : unit_factor = ref.unit_Velocity
		elif(field == 'vx3.') : unit_factor = ref.unit_Velocity
		elif(field == 'Tgas.') : unit_factor = ref.unit_Temperature
		elif(field == 'ionx.') : unit_factor = 1.0
		elif(field == 'iony.') : unit_factor = 1.0
		elif(field == 'prs.') : unit_factor = ref.unit_Pressure
		elif(field == 'ueuv.') : unit_factor = 1.0/ref.unit_Volume
		elif(field == 'urec.') : unit_factor = 1.0/ref.unit_Volume
		elif(field == 'fun.') : unit_factor = 1.0
		elif(field == 'bx1.') : unit_factor = ref.unit_MagneticField
		elif(field == 'bx2.') : unit_factor = ref.unit_MagneticField
		elif(field == 'bx3.') : unit_factor = ref.unit_MagneticField
		elif(field == 'phi.') : unit_factor = ref.unit_GravPotential
		else : raise ValueError('The field flag does not exist')

	data = data*unit_factor

	if file.tell() != eof: print ('Error: Too few bytes read.')

	file.close()
	return data

def readtimerange(filedirectory,field,N=200,tbeg=100,tend=100) :
	"""
		filedirectory: string 
			directory of files
		field: string
			The field to read
		N : integer
			Resolution. 200^3 by default
		tbeg: integer
			First timestep to read. 100 by default. 
		tend: integer 
			Last timestep to read. 100 by default. 
	"""

	print("Reading {} from directory {} for times {} to {}".format(field,filedirectory,tbeg,tend))
	i=tbeg
	data = np.zeros((tend-tbeg+1,N*N*N))
	while i<=tend:
		data[i-tbeg] = readsinglefile(filedirectory,timestep=tbeg+i,field=field,N=N)
		i = i+1

	print("Reading done")
	return data





