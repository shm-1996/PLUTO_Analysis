from header import *

def plotColumnDensity(directory,tstart=100,tend=300,N=200,outdir=None):

	directory = os.path.abspath(directory)
	m = constant.mH_HydrogenMass
	k_boltzmann=constant.k_BoltzmannConstant
	dl = (4.0/N)*ref.unit_Length
	time = tstart
	columndens_z = np.zeros((tend-tstart+1,N,N))
	columndens_y = np.zeros((tend-tstart+1,N,N))
	columndens_x = np.zeros((tend-tstart+1,N,N))

	#Compute Column Densities
	while time <=tend :
		rho1d    = read.readsinglefile(directory,time,N,'rho')*ref.unit_Density
		iongas1d = read.readsinglefile(directory,time,N,'ionx')
		mu1d = (iongas1d*0.5+(1.-iongas1d)*1.0)*m
		no_rho = np.reshape(rho1d/mu1d,(N,N,N))
		columndens_z[time-tstart] = no_rho.sum(0) *dl 
		columndens_x[time-tstart] = no_rho.sum(2) *dl 
		time = time+1

	max_columndens_z = np.max(np.log10(columndens_z))
	min_columndens_z = np.min(np.log10(columndens_z))
	max_columndens_x = np.max(np.log10(columndens_z))
	min_columndens_x = np.min(np.log10(columndens_z))

	if(outdir is None):
		outdir = directory + '/Column_Density/'
	else:
		outdir = os.path.abspath(outdir)

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	time = tstart


	while time <=tend :
		fig,axs = plt.subplots(ncols= 2,figsize=(14,12))
		fig.subplots_adjust(wspace=0.02)
		im = axs[0].imshow(np.log10(columndens_z[time-tstart]),extent=[0,4,0,4],cmap='Spectral_r',origin='lower',vmin=min_columndens_z,vmax=max_columndens_z)
		axs[0].set_xticks([0.5,1.5,2.5,3.5])
		axs[0].set_yticks([0.5,1.5,2.5,3.5])
		axs[0].tick_params(which='major',width=1.5,length=4)
		axs[0].tick_params(which='minor',width=0.7,length=2)
		axs[0].set_xlabel(r'$x\;$(pc)',labelpad=0.03)
		axs[0].set_ylabel(r'$y\;$(pc)',labelpad=0.03)
		time_sim = (time-tstart)*10
		axs[0].annotate(r"$\mathbf{t} = \bf{%d}\;\mathrm{\mathbf{Kyr}}$"%time_sim,(0.2,3.7),color='#3E403D',fontsize=12,weight='black',alpha=0.8)
		#cax = fig.add_axes([0.12,0.02,0.02,0.2])
		cbar = fig.colorbar(im,ax = axs[0],shrink=1.0,use_gridspec=True,orientation='horizontal',pad=0.05)
		cbar.outline.set_linewidth(2)
		cbar.ax.tick_params(which='major',width=1.5,length=4)
		cbar.ax.set_xlabel(r"$\log_{10}\mathrm{N}_z\;(\mathrm{cm}^{-2})$",rotation=0,labelpad=5,fontsize=16)

		im = axs[1].imshow(np.transpose(np.log10(columndens_x[time-tstart])),extent=[0,4,0,4],cmap='Spectral_r',origin='lower',vmin=min_columndens_x,vmax=max_columndens_x)
		axs[1].set_xticks([0.5,1.5,2.5,3.5])
		axs[1].set_yticks([0.5,1.5,2.5,3.5])
		axs[1].set_yticklabels([])
		axs[1].tick_params(which='major',width=1.5,length=4)
		axs[1].tick_params(which='minor',width=0.7,length=2)
		axs[1].set_xlabel(r'$z\;$(pc)',labelpad=0.03)
		axs[1].set_xlabel(r'$z\;$(pc)')
		time_sim = (time-tstart)*10
		axs[1].annotate(r"$\mathbf{t} = \bf{%d}\;\mathrm{\mathbf{Kyr}}$"%time_sim,(0.2,3.7),color='#3E403D',fontsize=12,weight='black',alpha=0.8)

		cbar = fig.colorbar(im,ax = axs[1],shrink=1.0,use_gridspec=True,orientation='horizontal',pad=0.05)
		cbar.outline.set_linewidth(2)
		cbar.ax.tick_params(which='major',width=1.5,length=4)
		cbar.ax.set_xlabel(r"$\log_{10}\mathrm{N}_x\;(\mathrm{cm}^{-2})$",rotation=0,labelpad=5,fontsize=16)
		plt.setp(axs[0].spines.values(),linewidth=2.0)
		plt.setp(axs[1].spines.values(),linewidth=2.0)
		plt.savefig(outdir+'%04d'%time+"CombinedColumn.eps",bbox_inches='tight')
		plt.close(fig)
		time = time+1

	print('Combined Column Densities Created \n')

	return

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
	args = vars(ap.parse_args())


	plotColumnDensity(args['directory'],args['tstart'],args['tend'],args['N'],args['outdir'])




