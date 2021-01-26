from header import *
import timeit
from joblib import Parallel, delayed
from itertools import repeat


#Some globals
m = constant.mH_HydrogenMass
k_boltzmann=constant.k_BoltzmannConstant
use_Density_Threshold = True
dens_threshold = 2.0
mass_weighted = True
fit_lognormal_M = False
fit_Hopkins = False


def hk(s,sigma_s,theta) :                # Function from Hopkins(2013b): <Mass Weighted Version     
    #sigma_s = -2*np.mean(s)*(1.+theta)
    lamb = (sigma_s**2)/(2*(theta**2))
    omega = lamb/(1+theta) - s/theta
    return (iv(1,2*np.sqrt(lamb*omega))*np.exp(-(lamb+omega))*np.sqrt(lamb/((theta**2) * omega)))*np.exp(s)

def b_method1(sigma_s,M) :
    return np.sqrt(np.exp(sigma_s**2)-1.0)/np.sqrt(M**2)

def std_error(sigma_s,M,error_s,error_M) :
	db_ds = sigma_s * np.exp(sigma_s**2)/(M*np.sqrt(np.exp(sigma_s**2) - 1.0))    
	db_dM = np.sqrt(np.exp(sigma_s**2) - 1.0)/(M**2)
	error_b = np.sqrt((db_ds*error_s)**2 + (db_dM*error_M)**2)
	return error_b

def gaus(x,a,x0,sigma):
		    return a*np.exp(-(x-x0)**2/(2*sigma**2))*np.exp(x)	

def calculate_b_direct(directory,time,N,Use_temp=False) : 
	rho1d	= read.readsinglefile(directory,time,N,'rho')
	if(Use_temp is True) :
		Tgas1d    = read.readsinglefile(directory,time,N,'Tgas')
		Tgas1d = np.mean(Tgas1d)
	else :
		Tgas1d = 10.
	vx1d    = read.readsinglefile(directory,time,N,'vx1')
	vy1d    = read.readsinglefile(directory,time,N,'vx2')
	vz1d    = read.readsinglefile(directory,time,N,'vx3')

	vcom_x = np.average(vx1d,weights = rho1d)
	vcom_y = np.average(vy1d,weights = rho1d)
	vcom_z = np.average(vz1d,weights = rho1d)
	vx = vx1d-vcom_x
	vy = vy1d-vcom_y
	vz = vz1d-vcom_z
	vel = np.sqrt(vx**2 + vy**2 + vz**2)
	cs = np.sqrt(k_boltzmann*Tgas1d/m)

	Mach_No = np.sqrt(np.mean(vel**2))/np.mean(cs)
	scale_dens = np.log(rho1d/np.mean(rho1d))
	sigma_scaled = np.std(rho1d/np.mean(rho1d))

	b = sigma_scaled/Mach_No
	return sigma_scaled, Mach_No, b

def calculate_b(directory,time,N,plot=False,outdir=None) :
	#Read data
	tbeg = time
	tend = tbeg 
	rho1d	= read.readsinglefile(directory,time,N,'rho')
	iongas1d = read.readsinglefile(directory,time,N,'ionx')
	Tgas1d    = read.readsinglefile(directory,time,N,'Tgas')
	vx1d    = read.readsinglefile(directory,time,N,'vx1')
	vy1d    = read.readsinglefile(directory,time,N,'vx2')
	vz1d    = read.readsinglefile(directory,time,N,'vx3')

	
	#Choose only neutral gas
	indices = np.where(iongas1d<1.e-7)

	#Compute COM velocity to subtract
	vcom_x = np.average(vx1d[indices[0]],weights = rho1d[indices[0]])
	vcom_y = np.average(vy1d[indices[0]],weights = rho1d[indices[0]])
	vcom_z = np.average(vz1d[indices[0]],weights = rho1d[indices[0]])
	vx_neutral = vx1d[indices[0]]-vcom_x
	vy_neutral = vy1d[indices[0]]-vcom_y
	vz_neutral = vz1d[indices[0]]-vcom_z


	#Compute some more quantities
	mu1d = (iongas1d[indices[0]]*0.5+(1.-iongas1d)[indices[0]]*1.0)*m
	cs_neutral = np.sqrt(k_boltzmann*Tgas1d[indices[0]]/mu1d)
	vel_neutral = np.sqrt(vx_neutral**2+vy_neutral**2+vz_neutral**2)
	rho_neutral = rho1d[indices[0]]
	nodens_neutral = np.log10(rho_neutral/(constant.mH_HydrogenMass))
	rho_mean = np.mean(rho_neutral)


	#Prepare scatter plot variables : s and log_10(Mach)
	s = np.log(rho_neutral/rho_mean)
	Mach_No = vel_neutral/cs_neutral
	if(mass_weighted) :
	    hist_2D,xedges,yedges = np.histogram2d(s, Mach_No, (60, 60),normed=True,weights=rho_neutral)
	else :
	    hist_2D,xedges,yedges = np.histogram2d(s, Mach_No, (60, 60),normed=True)
	xidx = np.clip(np.digitize(s, xedges), 0, hist_2D.shape[0]-1)
	yidx = np.clip(np.digitize(Mach_No, yedges), 0, hist_2D.shape[1]-1)
	color = hist_2D[xidx, yidx]
	no_dens = np.log10(rho_neutral/(constant.mH_HydrogenMass))

	#Apply density threshold
	if(use_Density_Threshold) :
	    indices_density = np.where(no_dens>dens_threshold)
	else :
	    threshold = 0.5*np.max(color)
	    indices_density = np.where(color>threshold)


	#Prepare bins
	density_bins = np.arange(np.floor(np.min(s)),np.ceil(np.max(s)),0.5)
	density_bins_threshold = np.arange(np.floor(np.min(s[indices_density[0]])),
		np.ceil(np.max(s[indices_density[0]])),0.5)

	Mach_bins = np.logspace(np.log10(np.min(Mach_No)),np.log10(np.max(Mach_No)),50)
	Mach_bins_threshold = np.logspace(np.log10(np.min(Mach_No[indices_density[0]])),
		np.log10(np.max(Mach_No[indices_density[0]])),50)


	if(fit_Hopkins) :

		RHO_DENSEMEAN = np.mean(rho_neutral[indices_density[0]])
		dense_gas = np.log(rho_neutral[indices_density[0]]/RHO_DENSEMEAN)
		hist = np.histogram(dense_gas,bins=density_bins_threshold,
			normed=True,weights=rho_neutral[indices_density[0]])
		b = hist[1]
		sigma_s = np.std(dense_gas)
		bins = (b[1:] + b[:-1])/2
		try :

			popt , pcov = curve_fit(hk,bins,hist[0],p0=(sigma_s,0.1))
			sigma_s = popt[0]
			T = popt[1]
			perr = np.sqrt(np.diag(pcov))
			error_sigma = perr[0]
		except :
		
			sigma_s = 0.0
			error_sigma = 100.0
			print('Fitting Hopkins was not possible at t=%d'%((time-tstart)*10))	

	else :
		
		hist = np.histogram(s[indices_density[0]],bins=density_bins_threshold,
			density=True,weights=rho_neutral[indices_density[0]])
		b = hist[1]
		bins = (b[1:] + b[:-1])/2					   
		popt,pcov = curve_fit(gaus,bins,hist[0],
			p0=[0.5,np.median(s[indices_density[0]]),np.std(s[indices_density[0]])])

		scale_dens = popt[0]  
		dens_0 = popt[1]
		sigma_s = np.abs(popt[2])
		perr = np.sqrt(np.diag(pcov))
		error_sigma = perr[2]
			

	if(fit_lognormal_M) :
		hist = np.histogram(Mach_No[indices_density[0]],bins=Mach_bins_threshold,density=True,weights=rho_neutral[indices_density[0]])
		b = hist[1]
		bins = (b[1:] + b[:-1])/2
		popt,pcov = curve_fit(gaus,bins,hist[0],p0=[0.5,np.median(Mach_No[indices_density[0]]),np.std(Mach_No[indices_density[0]])])
		scale_Mach = popt[0]  
		Mach_0 = popt[1]
		sigma_Mach = np.abs(popt[2])
		perr = np.sqrt(np.diag(pcov))
		error_Mach = perr[2]

	else :
		Mach = Mach_No[indices_density[0]]
		Dens = rho_neutral[indices_density[0]]
		#Dens= Dens*0.0 + 1.0
		M_sqr = np.sum(Mach**2 * Dens)/np.sum(Dens)
		M_lin = np.sum(Mach * Dens)/np.sum(Dens)
		sigma_Mach = (M_sqr-M_lin**2)**0.5
		error_Mach = 0.0
		

	b1 = b_method1(sigma_s,sigma_Mach)

	if(plot is True):
		if(outdir is None):
			outdir = os.path.abspath(directory) + '/Analysis/PDF_Scatters/'
		else:
			outdir = os.path.abspath(outdir) + '/PDF_Scatters/'

		if not os.path.exists(outdir):
			os.makedirs(outdir)	


		plt.clf()
		plt.style.use('classic')


		# Sum of counts of Total Histogram : s 
		normalisation_constant_x = np.sum(np.histogram(s,bins=density_bins,
			weights=rho_neutral)[0]*np.diff(np.histogram(s,bins=density_bins,weights=rho_neutral)[1]))
		
		# Sum of counts of Dense Gas Histogram : s
		normalisation_dense_x = np.sum(np.histogram(s[indices_density[0]],
			bins=density_bins_threshold,
			weights=rho_neutral[indices_density[0]])[0]*np.diff(np.histogram(s[indices_density[0]],
				bins=density_bins_threshold,weights=rho_neutral[indices_density[0]])[1]))

		# Sum of counts of Total Histogram : Mach 
		normalisation_constant_y = np.sum(np.histogram(Mach_No,bins=Mach_bins,
			weights=rho_neutral)[0]*np.diff(np.histogram(Mach_No,bins=Mach_bins,
				weights=rho_neutral)[1]))


		# Sum of counts of Dense Gas Histogram : Mach
		normalisation_dense_y = np.sum(np.histogram(Mach_No[indices_density[0]],
			bins=Mach_bins_threshold,
			weights=rho_neutral[indices_density[0]])[0] * np.diff(np.histogram(Mach_No[indices_density[0]],
				bins=Mach_bins_threshold,weights=rho_neutral[indices_density[0]])[1]))

		#Start Plot
		fig, axs = plt.subplots(2,2,gridspec_kw={'height_ratios': [1,2],'width_ratios':[2,1.5]},figsize=(10,10))
		fig.subplots_adjust(hspace=0.0)
		fig.subplots_adjust(wspace=0.0)


		axs[1,0].set_xlim(-14.8,8.)
		axs[1,0].set_ylim(1.e-2,8.e2)
		im = axs[1,0].scatter(s, Mach_No, c=color, s=20, alpha=0.5, edgecolors=None,cmap='rainbow')
		cax = fig.add_axes([0.2,0.15,0.02,0.15])
		cbar = fig.colorbar(im,ax = axs[1],cax=cax,orientation='vertical',shrink=0.5,ticks=[0.01,0.02,0.03])
		if(mass_weighted) :
		    cax.set_ylabel(r"$\bf{P_{\mathrm{MW}}}$",rotation=0,fontsize=14,labelpad=-70)
		else :
		    cax.set_ylabel(r"$\bf{P_V}$",rotation=0,fontsize=14,labelpad=-60)

		axs[1,0].tick_params(axis='x')
		axs[1,0].set_ylabel(r'Local Mach number ($M$)',fontsize=16)
		axs[1,0].set_xlabel(r'$s = \ln(\rho/\rho_0$)',fontsize=18)
		axs[1,0].tick_params(which='major',width=1.5,length=4)
		axs[1,0].tick_params(which='minor',width=0.7,length=2)
		if(use_Density_Threshold) :
		    axs[1,0].axvline(np.min(s[indices_density[0]]),color='r',linestyle='--',linewidth = 1, label='Density Threshold')
		else :
		    axs[1,0].tricontour(s, Mach_No, color,(0.0,threshold),alpha=1.0,colors='k',linewidths=1.5)
		    
		axs[0,0].get_shared_x_axes().join(axs[1,0], axs[0,0])
		hist,bins = np.histogram(s[indices_density[0]],bins=density_bins_threshold,
			weights=rho_neutral[indices_density[0]])
		widths = np.diff(bins)[0]
		axs[0,0].set_ylabel(r'$\log_{10}\;P_{\mathrm{MW}}(s)$',fontsize=18)
		plot_bins = (bins[1:] + bins[:-1])/2.
		axs[0,0].bar(plot_bins, hist/normalisation_constant_x, widths,log=True,edgecolor='k', color='skyblue',lw=1.5,fill=True)
		axs[0,0].hist(s,bins=density_bins,density=True,log=True, edgecolor = 'black',weights=rho_neutral,lw=0.1,histtype='stepfilled',color='#76748090')
		#axs[0,0].plot(bins,hk(bins,popt[0])*(normalisation_dense_x/normalisation_constant_x),'.-',label=r'fit: $\sigma_s$ = %0.2f $\pm$ %0.2f'%(popt[0],perr))
		axs[0,0].plot(plot_bins,gaus(plot_bins,scale_dens,dens_0,sigma_s)*(normalisation_dense_x/normalisation_constant_x),'-',lw=4.0,color='#FF5D00',label=r'fit: $\sigma_s$ = %0.2f $\pm$ %0.2f'%(sigma_s,error_sigma))
		axs[0,0].plot(plot_bins,gaus(plot_bins,scale_dens,dens_0,sigma_s)*(normalisation_dense_x/normalisation_constant_x),'-',lw=1.0,color='w')
		axs[0,0].legend(frameon=False,loc='upper left',fontsize=15)
		axs[0,0].set_xlim(-14.8,8.)
		limits = axs[0,0].get_ylim()
		axs[0,0].set_ylim(5.e-5,limits[1])
		axs[0,0].set_xticklabels([])
		axs[0,0].set_xticks([])
		axs[0,0].tick_params(which='major',width=1.5,length=4)
		axs[0,0].tick_params(which='minor',width=0.7,length=2)


		ax2 = axs[0,0].twiny()  # instantiate a second axes that shares the same x-axis
		ax2.set_xlabel(r'$\log_{10}\; n \;(\mathrm{cm}^{-3})$',fontsize=18)  # we already handled the x-label with ax1
		#ax2.scatter(no_dens,Mach_No,c=color,s=0,edgecolor='',cmap='rainbow',alpha=0.0)
		s_no_convert = lambda s : np.log10(np.exp(s) *(rho_mean/(m)))
		ax2.tick_params(which='major',width=1.5,length=4)
		s_min,s_max = axs[0,0].get_xlim()
		ax2.set_xlim((s_no_convert(s_min)),(s_no_convert(s_max)))
		ax2.plot([],[])
		ax2.set_xticks((-2,0,2,4,6))
		ax2.tick_params(axis='x') 
		if(use_Density_Threshold) :
		    axs[0,0].axvline(np.min(s[indices_density[0]]), color='r', linestyle='--', linewidth=1,label='Threshold Density MC')



		axs[1,1].get_shared_y_axes().join(axs[1,0], axs[1,1])

		hist,bins = np.histogram(Mach_No[indices_density[0]],bins=Mach_bins_threshold,weights=rho_neutral[indices_density[0]])
		widths = np.diff(bins)

		axs[1,1].set_xlabel(r'$\log_{10}\; P_\mathrm{MW}(M)$',fontsize=18)
		axs[1,1].tick_params(which='major',width=1.5,length=4)
		axs[1,1].tick_params(which='minor',width=0.7,length=2)

		axs[1,1].barh(np.sqrt(bins[1:] * bins[:-1]), hist/(normalisation_constant_y*widths), widths,log=True,edgecolor='k', color='skyblue',lw=1.5,fill=True)
		axs[1,1].hist(Mach_No,bins=Mach_bins,density=True,log=True,histtype='stepfilled',color='#76748090',orientation='horizontal',edgecolor='black',weights=rho_neutral,lw=0.1)
		if(fit_lognormal_M) :
		    axs[1,1].plot(gaus(bins,scale_Mach,Mach_0,sigma_Mach),bins,'.-',label=r'fit: $\mathcal{M}$ = %0.2f $\pm$ %0.2f'%(sigma_Mach,error_Mach),color='k',linewidth=2)
		else :
		    axs[1,1].text(0.06,0.05,r'$\mathcal{M} = \langle M_{\mathrm{MW}}^2 \rangle - {\langle M_{\mathrm{MW}}\rangle }^2 = %0.2f$'%sigma_Mach,transform=axs[1,1].transAxes,fontsize=16)
		limits = axs[1,1].get_xlim()
		axs[1,1].set_xlim(limits[0]*5000,limits[1])
		#axs[1,1].set_ylim([limits_0[0],limits_0[1]+0.3])
		#axs[1,1].axhline(MPIL_MEDIAN, color='r', linestyle='--', linewidth=1,label=r'$\log_{10}(\mathcal{M}) = %0.2f \pm %0.2f$'%(MPIL_MEDIAN,MPIL_STD))
		if(fit_lognormal_M):
			axs[1,1].legend(fontsize=12,loc=8,frameon=False)
		axs[1,1].set_yscale("log")
		axs[1,1].set_yticklabels([])
		axs[0,1].axis('off')
		axs[1,1].set_ylim(1.e-2,8.e2)
		plt.setp(axs[0,0].spines.values(),linewidth=2.0)
		plt.setp(axs[1,0].spines.values(),linewidth=2.0)
		plt.setp(axs[1,1].spines.values(),linewidth=2.0)
		plt.savefig(outdir+'Scatter_{:04d}.png'.format(time),dpi=500,bbox_inches='tight')
		#plt.savefig(outdir+'Scatter_{:04d}.pdf'.format(time),bbox_inches='tight')

	if(fit_lognormal_M) :
		return sigma_s,error_sigma,sigma_Mach,error_Mach,b1	
	else :
		return sigma_s,error_sigma,sigma_Mach,b1

def calculate_b_range(directory,N,tstart,tfinish,plot=False,outdir=None):
	print("Calculating b in time range {} to {}".format(tstart,tfinish))
	SIGMA_S = np.zeros(tfinish-tstart+1)
	ERROR_S = np.zeros(tfinish-tstart+1)
	MACH_NO = np.zeros(tfinish-tstart+1)
	if(fit_lognormal_M) :
		ERROR_M = np.zeros(tfinish-tstart+1)
	B_VALUE = np.zeros(tfinish-tstart+1)
	B_ERROR = np.zeros(tfinish-tstart+1)

	time = tstart
	#Parallelise this calculation
	time_arr = np.arange(tstart,tfinish+1,1)
	results = Parallel(n_jobs=-1,prefer='processes',verbose=1)(map(delayed(calculate_b),
		repeat(directory),time_arr,repeat(N),repeat(plot),repeat(outdir)))

	while time <=tfinish :
		if(fit_lognormal_M) :
			sigma_s,error_sigma,Mach_No,error_Mach,b1 = results[time-tstart]
		else :
			sigma_s,error_sigma,Mach_No,b1 = results[time-tstart]				
			error_Mach = 0.0
		SIGMA_S[time-tstart] = sigma_s
		ERROR_S[time-tstart] = error_sigma
		MACH_NO[time-tstart] = Mach_No
		if(fit_lognormal_M) :
			ERROR_M[time-tstart] = error_Mach
		B_VALUE[time-tstart] = b1 
		if(error_sigma==100.0) :
			B_ERROR = 100.0
		else :	
			B_ERROR[time-tstart] = std_error(sigma_s,Mach_No,error_sigma,error_Mach)
		time = time + 1

	if(fit_lognormal_M):
		return SIGMA_S,ERROR_S,MACH_NO,ERROR_M,B_VALUE,B_ERROR
	else:
		return SIGMA_S,ERROR_S,MACH_NO,B_VALUE,B_ERROR


def calculate_b_directrange(directory,N,tstart,tfinish):
	print("Calculating b in time range {} to {}".format(tstart,tfinish))
	SIGMA_S = np.zeros(tfinish-tstart+1)
	MACH_NO = np.zeros(tfinish-tstart+1)
	B_VALUE = np.zeros(tfinish-tstart+1)

	time = tstart
	while time <=tfinish :

		sigma_s,M,b1 = calculate_b_direct(directory,time,N)

		SIGMA_S[time-tstart] = sigma_s
		MACH_NO[time-tstart] = M 
		B_VALUE[time-tstart] = b1 
		time = time + 1
	return SIGMA_S,MACH_NO,B_VALUE




def plot_quantities(tfinish,Sigma_s,Mach_No,b,save=False,outdir=None):
	print("Plotting summarised b quantities.")

	fig,axs = plt.subplots(nrows=3,figsize=(8,6),sharex=True)
	time = np.linspace(0,tfinish,tfinish+1)*10.0
	axs[0].plot(time,Sigma_s,'.-')
	axs[0].set_ylabel(r'$\sigma_s$')

	axs[1].plot(time,Mach_No,'.-')
	axs[1].set_ylabel(r'$\mathcal{M}$')

	axs[2].plot(time,b,'.-')
	axs[2].set_xlabel(r'$ t \, (\mathrm{kyr})$')
	axs[2].set_ylabel(r'$b$')

	if(outdir is None):
		outdir = os.getcwd() + '/'
	if(save):
		plt.savefig(outdir+'b_summary',bbox_inches='tight')
	else:
		plt.show()

if __name__ == "__main__":

	ap = argparse.ArgumentParser(description='Command Line Inputs for b computing script.')
	ap.add_argument('-directory', metavar='Output Directory', action='store', type=str,
		default=os.path.realpath('./')+'/',help='Output files directory')
	ap.add_argument('-N', metavar='Resolution', action='store', type=int,default=200,help='Resolution of Simulation',dest='N')
	ap.add_argument('-output', metavar='Output Directory of table', action='store', type=str,
	                default=None,help='Directory to save table.')
	ap.add_argument('-tstart', metavar='t0', action='store', type=int,default=100,
	            help='Starting timestep.')
	ap.add_argument('-tend', metavar='tend', action='store', type=int,default=300,
	            help='Last timestep.')
	ap.add_argument('-direct', action='store_true',
					help='Calculate b from direct formulation without fitting assuming isothermal.')
	ap.add_argument('-table',action='store_true',
				help='Flag to write b to table. Default false.')
	ap.add_argument('-plot',action='store_true',
		help='Flag to plot scatter plot of s and Mach Number for each timestep.')
	args = vars(ap.parse_args())

	# time the script
	start_time = timeit.default_timer()

	directory = os.path.abspath(args['directory']) + '/'
	if(args['output'] is None):
		output = directory + 'Analysis/'
	else:
		output = os.path.abspath(args['output']) + '/'
	N = args['N']
	tstart = args['tstart']
	tfinish = args['tend']
	
	if not os.path.exists(output):
		os.makedirs(output)	

	if(args['direct'] is True):
		Sigma_s, Mach_No, b = calculate_b_directrange(directory,N,tstart,tfinish)
		pkl_obj = [Sigma_s, Error_s, Mach_No, b]
	else:
		Sigma_s, Error_s, Mach_No, b, Error_b = calculate_b_range(directory,N,tstart,tfinish,
			plot=args['plot'],outdir=output)
		pkl_obj = [Sigma_s, Error_s, Mach_No, b, Error_b]

	saveObj(pkl_obj,output+'b_summary')
	
	if(args['table'] is True):
		print("Writing to table")
		# Write Table of b
		time_column = np.arange((tstart-100)*10,(tfinish-100+1)*10,10)
		if(fit_lognormal_M) :
			header = '1. Time(Kyr) \t sigma_s \t Error_s \t Mach No \t Error_M \t b \t error_b \n\n'
			text_columns = np.column_stack((time_column,Sigma_s,Error_s,Mach_No,ERROR_M,b,Error_b))
			np.savetxt(output+'b_timevariation_LognormalM.dat',text_columns,fmt='%0.2f',delimiter='\t',header=header)
		elif(args['direct'] is True) :
			header = '1. Time(Kyr) \t sigma_s \t Error_s \t Mach No \t b \t error_b \n\n'
			text_columns = np.column_stack((time_column,Sigma_s,Mach_No,b))
			np.savetxt(output+'b_direct.dat',text_columns,fmt='%0.2f',delimiter='\t',header=header)
		else:
			header = '1. Time(Kyr) \t sigma_s \t Error_s \t Mach No \t b \t error_b \n\n'
			text_columns = np.column_stack((time_column,Sigma_s,Error_s,Mach_No,b,Error_b))
			np.savetxt(output+'b_timevariation.dat',text_columns,fmt='%0.2f',delimiter='\t',header=header)	
	
	plot_quantities(tfinish-tstart,Sigma_s,Mach_No,b,save=True,outdir=output)

	print('Time Variation of b Calculated\n\n')	

	# time the script
	stop_time = timeit.default_timer()
	total_time = stop_time - start_time
	print("***************** time to finish = "+str(total_time)+"s *****************")




