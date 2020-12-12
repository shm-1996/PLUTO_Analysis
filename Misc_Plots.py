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
        axs.set_ylabel(r'$p_V \, \log_{10} \, \left( \rho \right)$')
        axs.set_xlabel(r'$\log_{10} \, \left( \rho \right)$')
        plt.savefig(outdir+'Density_%04d'%time,bbox_inches='tight')
        plt.clf()
        plt.close()

    return



    

    

