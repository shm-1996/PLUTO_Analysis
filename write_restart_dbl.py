from numpy import *
import numpy as np 
import os
import argparse
import sys
from shutil import copyfile
import glob

ap = argparse.ArgumentParser(description='Command Line Inputs for restartings. All inputs optional. ')
ap.add_argument('-directory','--directory',metavar='Output_Directory',default=None,
	help='an argument for specifying the output directory.',type=str)
ap.add_argument('-directory2','--directory2',metavar='Copy_directory',
	default=None,help='an argument for specifying the output directory to copy from.',type=str)

ap.add_argument('-N', metavar='Resolution', action='store', type=int,default=200,
                help='Resolution of Simulation. Default = 200',dest='N')
ap.add_argument('-t', metavar='Timestep', action='store', type=int,default=100,
                help='Timestep of Simulation Output. Default = 100',dest='time')
args = vars(ap.parse_args())

def create_files(directory,N=200,time=100) :

	ionx = np.full(N*N*N,1.e-40)
	file = open(directory+"ionx.{:04d}.dbl".format(time),"wb")
	file.write(ionx)
	file.close()

	ueuv = np.full(N*N*N,1.e-40)
	file = open(directory+"ueuv.{:04d}.dbl".format(time),"wb")
	file.write(ueuv)
	file.close()

	urec = np.full(N*N*N,1.e-40)
	file = open(directory+"urec.{:04d}.dbl".format(time),"wb")
	file.write(urec)
	file.close()

	phi = np.full(N*N*N,1.e-40)
	file = open(directory+"phi.{:04d}.dbl".format(time),"wb")
	file.write(phi)
	file.close()

	#Tgas 
	Tgas = np.full(N*N*N,10.0)
	file = open(directory+"Tgas.{:04d}.dbl".format(time),"wb")
	file.write(Tgas)
	file.close()

	#iony
	iony = np.full(N*N*N,1.0)
	file = open(directory+"iony.{:04d}.dbl".format(time),"wb")
	file.write(iony)
	file.close()


def copy_files(input_directory,output_directory,fileno):
	#Copying dbl files
	for inputfile in glob.glob(input_directory+'*.{:04d}.dbl'.format(fileno)):
		file_suffix = inputfile.split('/')[-1]
		copyfile(inputfile,output_directory+file_suffix)

	#Copying .out files

	for inputfile in glob.glob(input_directory+'*.out'):
		file_suffix = inputfile.split('/')[-1]
		copyfile(inputfile,output_directory+file_suffix)

	return

	return
if __name__ == '__main__':
	if(not args['directory']) :
		print("Output Directory Required")
		sys.exit()
	print("Creating restart appropriate dbl files for a {}^3 resolution simulation in {}".format(args['N'],args['directory']))
	directory = os.path.abspath(args['directory']) + '/'
	create_files(directory,args['N'],args['time'])

	if(args['directory2']):
		directory2 = os.path.abspath(args['directory2']) + '/'
		print("Copying files for restart from {} to {}".format(directory2,directory))
		copy_files(directory2,directory,args['time'])

