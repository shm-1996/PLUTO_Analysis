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
        data = read.readsinglefile(directory,time,N,field)

        fig,axs = plt.subplots(ncols=1)
        axs[0].hist(np.log10(dens),bins=20,density=True,color="#40D3FF30",histtype="bar")
        axs[0].set_xlim(np.log10(min_dens),np.log10(max_dens))
        axs[0].set_ylim(r'$p_V \, \log_{10} \, \left( \rho \right)$')
        axs[0].set_xlim(r'$\log_{10} \, \left( \rho \right)$')
        plt.savefig(outdir+'Density_%04d'%time,bbox_inches='tight')
        plt.clf()
        plt.close()

    return



    

    

