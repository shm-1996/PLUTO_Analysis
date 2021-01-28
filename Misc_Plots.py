from header import *
from FunctionInterface import *
from Power_Spectra import *

fields = np.array(['rho','vx1','vx2','vx3','prs','Tgas','ionx','iony','ueuv','urec',
            'fun','bx1','bx2','bx3','phi'])

labels = np.array([r'$\rho$',r'$v_{\mathrm{x}}$',r'$v_{\mathrm{y}}$',r'$v_{\mathrm{z}}$',
r'$P_{\mathrm{th}}$',r'$T$',r'$x_{\mathrm{ion}}$',r'$y_{\mathrm{ion}}$',r'$u_{\mathrm{EUV}}$',
r'$u_{\mathrm{REC}}$','fun',r'$B_{\mathrm{x}}$',r'$B_{\mathrm{y}}$',r'$B_{\mathrm{z}}$',
r'$\phi$'])

unit_string = np.array([r'$\mathrm{g}\,\mathrm{cm}^{-3}$',r'$\mathrm{km}\,\mathrm{s}^{-1}$',
    r'$\mathrm{km}\,\mathrm{s}^{-1}$',r'$\mathrm{km}\,\mathrm{s}^{-1}$',
r'$\mathrm{erg} \, \mathrm{cm}^{-3}$',r'$\mathrm{K}$',None,None,
r'$\mathrm{photons} \, \mathrm{cm}^{-3}$',r'$\mathrm{photons} \, \mathrm{cm}^{-3}$',
'',r'$\mu \mathrm{G}$',r'$\mu \mathrm{G}$',r'$\mu \mathrm{G}$',
r'$\mathrm{erg} \, \mathrm{cm}^{-3}$'])

unit_factor = np.array([1.0,1.e-5,1.e-5,1.e-5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    1.e-6,1.e-6,1.e-6])

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


def Slice_Plot(directory,field='rho',tstart=100,tend=300,N=200,outdir=None,log=True,
    slice_index=-1):

    """
    Routine to plot slice plots of a field in a time range
        directory : string
            directory where files reside
        field: string
            Data field to read in
        tstart: integer 
            directory in which files reside
        tend : integer 
            timestep of file        
        N : integer
            Resolution
        outdir : string 
            Location to store the files.
        log : boolean
            Flag to use logscale in plots. True by default. 
        slice_index : integer
            Optional argument to compute min max in a slice of the box
        
    """

    directory = os.path.abspath(directory)
    m = constant.mH_HydrogenMass
    k_boltzmann=constant.k_BoltzmannConstant
    time = tstart

    if(outdir is None):
        outdir = directory + '/Slice_Plots/'
    else:
        outdir = os.path.abspath(outdir) + '/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    min_data,max_data = Compute_min_max(tstart,tend,directory,N,field,
        slice_index=slice_index)
    if(log):
        min_data = np.log10(min_data)
        max_data = np.log10(max_data)

    index = np.where(fields==field)[0]
    label_plot = labels[index][0]
    unit_plot = unit_factor[index][0]
    if(slice_index == -1):
        slice_index = int(N/2)+1
    label_plot = "{}".format(labels[index][0])
    if(log):
        label_plot = r"$\log_{10} \,$" +"{}".format(labels[index][0])

    #Add unit string
    if(unit_string[index][0]):
        label_plot = label_plot + r'$\;$' + '({})'.format(unit_string[index][0])

    for time in range(tstart,tend+1):
        data = read.readsinglefile(directory,time,N,field).reshape(N,N,N)
        #Central slice
        data = data[int(N/2)+1]
        data = data*unit_plot
        fig,axs = plt.subplots(ncols=1)
        if(log):
            im = axs.imshow(np.log10(data),vmin=min_data,vmax=max_data,extent=[0,4,0,4],
                cmap='Spectral_r',origin='lower')
        else:
            im = axs.imshow(data,vmin=min_data,vmax=max_data,extent=[0,4,0,4],
                cmap='Spectral_r',origin='lower')
        axs.set_xticks([0.5,1.5,2.5,3.5])
        axs.set_yticks([0.5,1.5,2.5,3.5])
        axs.set_xlabel(r'$x\;$(pc)',labelpad=0.03)
        axs.set_ylabel(r'$y\;$(pc)',labelpad=0.03)
        time_sim = (time-tstart)*10
        axs.annotate(r"$\mathbf{t} = \bf{%d}\;\mathrm{\mathbf{Kyr}}$"%time_sim,(0.2,3.7),
            color='#3E403D',fontsize=12,weight='black',alpha=0.8)
        cbar = fig.colorbar(im,ax = axs,shrink=1.0,use_gridspec=True,
                        orientation='vertical',pad=0.02)
        cbar.ax.set_ylabel(label_plot,rotation=90,
                   labelpad=5,fontsize=16)
        plt.savefig(outdir+'{}_{:04d}'.format(field,time),bbox_inches='tight')
        plt.clf()
        plt.close(fig)

    return


def Initial_Analyses(directory,tstart=100,tend=300,N=200):

    print("Performing initial analyses for this simulation....")
    print("Density PDF..")
    Density_PDF(directory,tstart,tend,N)

    print("Column Density...")
    plotColumnDensity(directory,tstart,tend,N)

    print("Making Movie")
    Make_Movie(directory,'CombinedColumn',tstart,tend,convert_pdfs=True)

    print("Power spectra..")
    Compute_Power_Spectra(directory,tstart,tend,N)

    print("Plotting Power spectra..")
    Plot_PowerSpectra(directory,tstart, tend, N)