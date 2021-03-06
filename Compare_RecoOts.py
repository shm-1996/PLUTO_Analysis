from pluto_header import *
from PlutoInterface import *
from Pluto_Plots import *

seedDirs = ["seed_120963","seed_130696","seed_140281",
    "seed_270797","seed_290691","seed_60468"]

OTS_directory = '/scratch/ek9/sm5890/Pillar_Simulations/OTS_Seeds'
Reco_directory = '/scratch/ek9/sm5890/Pillar_Simulations/RecombinationField_Seeds'        

def CompareOTSReco(seed=120963,field='rho',tstart=0,tend=100,N=200,
    outdir=None,log=True,slice_index=-1,show=False,projection=False):
    """
    Routine to plot compare slice plots across seeds
            basedirectory : string
                base directory where subdirectories reside
            field: string
                Data field to read in
            tstart: integer 
                directory in which files reside
            tend : integer 
                timestep of file        
            outdir : string 
                Location to store the files.
            log : boolean
                Flag to use logscale in plots. True by default. 
            slice_index : integer
                Optional argument to compute min max in a slice of the box
            show : Boolean
                Optional argument to only show not save. Can only be used 
                if tstart = tend
            projection: Boolean
                Optional argument to choose projection instead of slice
    
    """

    if(tend>tstart and show is True):
        raise ValueError("The 'show' argument can only be used if tstart=tend.")

    basedirectory = '/scratch/ek9/sm5890/Pillar_Simulations/'
    if(outdir is None):
        outdir = basedirectory + '/Compare_RecoOts/'
    else:
        outdir = os.path.abspath(outdir) + '/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    min_data = np.inf
    max_data = -np.inf
    #Find minimum/maximum across the 6 runs
    Dirs = [Reco_directory + "/seed_{}".format(seed),
            OTS_directory + "/seed_{}".format(seed)]
    for seedDir in Dirs:
        directory = seedDir
        min_data_temp,max_data_temp = Compute_min_max(tstart,tend,directory,N,field,
            slice_index=slice_index)
        if(log):
            min_data_temp = np.log10(min_data_temp)
            max_data_temp = np.log10(max_data_temp)
        if(max_data_temp>max_data):
            max_data = max_data_temp
        if(min_data_temp<min_data):
            min_data = min_data_temp

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
        fig,axs = plt.subplots(ncols=2,figsize=(8,4))
        i = 0
        for seedDir in Dirs:
            directory = seedDir
            data = read.readsinglefile(directory,time,N,field).reshape(N,N,N)
            data = data[int(N/2)+1]
            data = data*unit_plot
            if(log):
                im = axs[i].imshow(np.log10(data),vmin=min_data,vmax=max_data,extent=[0,4,0,4],
                    cmap='Spectral_r',origin='lower')
            else:
                im = axs[i].imshow(data,vmin=min_data,vmax=max_data,extent=[0,4,0,4],
                    cmap='Spectral_r',origin='lower')
            axs[i].set_xticks([0.5,1.5,2.5,3.5])
            axs[i].set_yticks([0.5,1.5,2.5,3.5])
            axs[i].set_xlabel(r'$x\;$(pc)',labelpad=0.03)
            if(i == 0):
                axs[i].set_ylabel(r'$y\;$(pc)',labelpad=0.03)

            if(i==1):
                cbar = fig.colorbar(im,ax = axs[i],shrink=1.0,use_gridspec=True,
                        orientation='vertical',pad=0.02)
                cbar.ax.set_ylabel(label_plot,rotation=90,
                           labelpad=5,fontsize=16)
            i = i+1

        time_sim = (time-tstart)*10
        axs[0].annotate(r"$\mathbf{t} = \bf{%d}\;\mathrm{\mathbf{Kyr}}$"%time_sim,(0.2,3.7),
            color='#3E403D',fontsize=12,weight='black',alpha=0.8)

        if(show):
            plt.show()                    
        else:
            plt.savefig(outdir+'{}{}_{:04d}'.format(seed,
                field,time),bbox_inches='tight')
            plt.clf()
            plt.close(fig)








        