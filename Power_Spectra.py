from header import *
from FunctionInterface import *
from Misc_Plots import * 
import subprocess

def Compute_PowerSpectra(directory,tstart=100,tend=300,N=200,outdir=None,plot=True):
    """
        Compute power spectra for a range of timesteps by calling external C++
        routine powspectrum.cpp. 
    """

    directory = os.path.abspath(directory) + '/'
    time = tstart

    if(outdir is None):
        outdir = directory + 'Power_Spectra/'
    else:
        outdir = os.path.abspath(outdir) + '/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #Compile Power spectra code 
    subprocess.run(["mpic++ -o pow powspectrum.cpp -lm -lfftw3_mpi -lfftw3"],shell=True)

    for time in range(tstart,tend+1):
        subprocess.run(["./pow {} {} {} {} {}".format(directory,time,N,N,N)],shell=True)

    if(plot is True):
        Plot_PowerSpectra(directory,tstart,tend,N)



