from header import *
from FunctionInterface import *

def Density_PDF(directory,tstart=100,tend=300,N=200,outdir=None):
    directory = os.path.abspath(directory)
    m = constant.mH_HydrogenMass
    k_boltzmann=constant.k_BoltzmannConstant
    time = tstart

    if(outdir is None):
        outdir = directory + '/PDFs/'
    else:
        outdir = os.path.abspath(outdir) + '/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    min_dens,max_dens = Compute_min_max(tstart,tend,directory,N,'rho')

    for time in range(tstart,tend+1):
        dens = read.readsinglefile(directory,time,N,'rho')

        fig,axs = plt.subplots(ncols=1)
        axs.hist(np.log10(dens),bins=40,density=True,color="#40D3FF90",histtype="bar")
        axs.set_xlim(np.log10(min_dens),np.log10(max_dens))
        axs.set_ylim(0.0,1.0)
        axs.set_ylabel(r'$p_V \, \log_{10} \, \left( \rho \right)$')
        axs.set_xlabel(r'$\log_{10} \, \left( \rho \right)$')
        plt.title(r'$t = {}$'.format(time*ref.unit_Time/const.Kyr))
        plt.savefig(outdir+'Density_%04d'%time,bbox_inches='tight')
        plt.clf()
        plt.close()

    return


def Plot_PowerSpectra(directory,tstart=100,tend=300,N=200):
    directory = os.path.abspath(directory) + '/'
    outdir = directory + 'Power_Spectra/'

    for time in range(tstart,tend):
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


    
    

    

