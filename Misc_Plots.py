from header import *
from FunctionInterface import *
from Power_Spectra import *

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

def Initial_Analyses(directory,tstart=100,tend=300,N=200):

    print("Performing initial analyses for this simulation....")
    print("Density PDF..")
    Density_PDF(directory,tstart,tend,N)

    print("Column Density...")
    plotColumnDensity(directory,tstart,tend,N)

    print("Making Movie")
    Make_Movie(directory,tstart,tend,convert_pdfs=True)

    print("Power spectra..")
    Compute_Power_Spectra(directory,tstart,tend,N)

    print("Plotting Power spectra..")
    Plot_PowerSpectra(directory,tstart, tend, N)