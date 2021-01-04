from header import *
from FunctionInterface import *
from Misc_Plots import * 
import subprocess
import socket

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


    if(socket.gethostname().split('-')[0] == 'gadi'):
        #Specific to gadi
        subprocess.run(["mpiCC -o pow powspectrum.cpp -lm -lfftw3_mpi -lfftw3"],shell=True)     
    else:
        subprocess.run(["mpic++ -o pow powspectrum.cpp -lm -lfftw3_mpi -lfftw3"],shell=True)

    for time in range(tstart,tend+1):
        subprocess.run(["./pow {} {} {} {} {}".format(directory,time,N,N,N)],shell=True)




def Plot_PowerSpectra(directory,tstart=100,tend=300,N=200):
    directory = os.path.abspath(directory) + '/'
    outdir = directory + 'Power_Spectra/'

    for time in range(tstart,tend+1):
        file = np.loadtxt(outdir+'_spect_vels%04d.dat'%time,skiprows=5)
        k = file[:,1]
        P_k = file[:,15]
        P_long = file[:,11]
        P_trv = file[:,13]

        fig,axs = plt.subplots(ncols=1)
        axs.plot(k,P_k,'x-',lw=2.0,label='Total')
        axs.plot(k,P_trv,'x-',lw=2.0,label='Transverse')
        axs.plot(k,P_long,'x-',lw=2.0,label='Longitudinal')
        axs.set_xlabel(r'$k$')
        axs.set_ylabel(r'$P(k)$')
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.legend(loc='best')
        plt.savefig(outdir+'Vels_%04d'%time,bbox_inches='tight')
        plt.clf()
        plt.close()

        file = np.loadtxt(outdir+'_spect_rhoweightedv%04d.dat'%time,skiprows=5)
        k = file[:,1]
        P_k = file[:,15]
        P_long = file[:,11]
        P_trv = file[:,13]

        fig,axs = plt.subplots(ncols=1)
        axs.plot(k,P_k,'x-',lw=2.0,label='Total')
        axs.plot(k,P_trv,'x-',lw=2.0,label='Transverse')
        axs.plot(k,P_long,'x-',lw=2.0,label='Longitudinal')
        axs.set_xlabel(r'$k$')
        axs.set_ylabel(r'$P(k)$')
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.legend(loc='best')
        plt.savefig(outdir+'RhoV_%04d'%time,bbox_inches='tight')
        plt.clf()
        plt.close()


if __name__ == "__main__":

    #Parsing Arguments
    ############################################################################
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for column density Plots. ')
    ap.add_argument('-directory',action='store',type=str,default='./',
        help='Directory of output files. By default current working directory." ')
    ap.add_argument('-outdir',action='store',type=str,default=None,
        help = 'Output directory for files.')
    ap.add_argument('-tstart',action='store',type=int,default=100,
        help = 'Start timestep.')
    ap.add_argument('-tend',action='store',type=int,default=300,
        help = 'End timestep.')
    ap.add_argument('-N',action='store',type=int,default=200,
        help = 'Resolution')
    ap.add_argument('-plot', action='store_true',
                    help='Flag to plot power spectra.')
    args = vars(ap.parse_args())

    Compute_PowerSpectra(args['directory'],args['tstart'],args['tend'],args['N'])
    if(args['plot'] is True):
        print("Plotting...")
        Plot_PowerSpectra(args['directory'],args['tstart'],args['tend'],args['N'])