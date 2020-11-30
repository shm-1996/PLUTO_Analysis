import numpy as np
import readdbl_timerange as read
import ReferenceUnits as ref
import matplotlib.pyplot as plt
import PhysicalConstantsCGS as constant 
import matplotlib.patches as patches
import os
import matplotlib as mpl
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
import sys
import time
from scipy.optimize import OptimizeWarning
import warnings
import argparse

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

def calculate_b_direct(directory,time,N) : 
	rho1d	= read.readsinglefile(directory,time,N,'rho')*ref.unit_Density
	try :
		Tgas1d    = read.readsinglefile(directory,time,N,'Tgas')
		Tgas1d = np.mean(Tgas1d)
	except :
		print("Setting Tgas to 10 K")
		Tgas1d = 10.
	vx1d    = read.readsinglefile(directory,time,N,'velx')*ref.unit_Velocity
	vy1d    = read.readsinglefile(directory,time,N,'vely')*ref.unit_Velocity
	vz1d    = read.readsinglefile(directory,time,N,'velz')*ref.unit_Velocity

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

def calculate_b(directory,time,N) :
	tbeg = time
	tend = tbeg 
	rho1d	= read.readsinglefile(directory,time,N,'rho')*ref.unit_Density
	iongas1d = read.readsinglefile(directory,time,N,'ionx')
	Tgas1d    = read.readsinglefile(directory,time,N,'Tgas')
	vx1d    = read.readsinglefile(directory,time,N,'velx')*ref.unit_Velocity
	vy1d    = read.readsinglefile(directory,time,N,'vely')*ref.unit_Velocity
	vz1d    = read.readsinglefile(directory,time,N,'velz')*ref.unit_Velocity

	

	indices = np.where(iongas1d<1.e-7)
	vcom_x = np.average(vx1d[indices[0]],weights = rho1d[indices[0]])
	vcom_y = np.average(vy1d[indices[0]],weights = rho1d[indices[0]])
	vcom_z = np.average(vz1d[indices[0]],weights = rho1d[indices[0]])
	vx_neutral = vx1d[indices[0]]-vcom_x
	vy_neutral = vy1d[indices[0]]-vcom_y
	vz_neutral = vz1d[indices[0]]-vcom_z


	# In[68]:


	mu1d = (iongas1d[indices[0]]*0.5+(1.-iongas1d)[indices[0]]*1.0)*m
	cs_neutral = np.sqrt(k_boltzmann*Tgas1d[indices[0]]/mu1d)
	vel_neutral = np.sqrt(vx_neutral**2+vy_neutral**2+vz_neutral**2)
	rho_neutral = rho1d[indices[0]]
	nodens_neutral = np.log10(rho_neutral/(constant.mH_HydrogenMass))
	RHO_MEAN = np.mean(rho_neutral)

	


	# In[72]:


	x = np.log(rho_neutral/RHO_MEAN)
	y = np.log10(vel_neutral/cs_neutral)
	Mach_No = vel_neutral/cs_neutral
	if(mass_weighted) :
	    hist_2D,xedges,yedges = np.histogram2d(x, y, (60, 60),normed=True,weights=rho_neutral)
	else :
	    hist_2D,xedges,yedges = np.histogram2d(x, y, (60, 60),normed=True)
	xidx = np.clip(np.digitize(x, xedges), 0, hist_2D.shape[0]-1)
	yidx = np.clip(np.digitize(y, yedges), 0, hist_2D.shape[1]-1)
	color = hist_2D[xidx, yidx]
	#idx = color.argsort()
	#x, y, color = x[idx], y[idx], color[idx]
	no_dens = np.log10(rho_neutral/(constant.mH_HydrogenMass))


	# In[73]:


	#indices_density = np.where((color>threshold) & (no_dens>2.))
	threshold = 0.5*np.max(color)
	
	if(use_Density_Threshold) :
	    indices_density = np.where(no_dens>dens_threshold)
	else :
	    indices_density = np.where(color>threshold)


	# In[74]:

	density_bins = np.arange(np.floor(np.min(x)),np.ceil(np.max(x)),0.1)
	density_bins_threshold = np.arange(np.floor(np.min(x[indices_density[0]])),np.ceil(np.max(x[indices_density[0]])),0.1)

	Mach_bins = np.arange(np.floor(np.min(y)),np.ceil(np.max(y)),0.1)
	Mach_bins_threshold = np.arange(np.floor(np.min(y)),np.ceil(np.max(y)),0.1)


	#plt.clf()

	if(fit_Hopkins) :

		RHO_DENSEMEAN = np.mean(rho_neutral[indices_density[0]])
		dense_gas = np.log(rho_neutral[indices_density[0]]/RHO_DENSEMEAN)
		hist = np.histogram(dense_gas,bins=density_bins_threshold,normed=True,weights=rho_neutral[indices_density[0]])
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
		
		hist = np.histogram(x[indices_density[0]],bins=density_bins_threshold,density=True,weights=rho_neutral[indices_density[0]])
		b = hist[1]
		bins = (b[1:] + b[:-1])/2

		
			   
		popt,pcov = curve_fit(gaus,bins,hist[0],p0=[0.5,np.median(x[indices_density[0]]),np.std(x[indices_density[0]])])

		scale_dens = popt[0]  
		dens_0 = popt[1]
		sigma_s = np.abs(popt[2])
		perr = np.sqrt(np.diag(pcov))
		error_sigma = perr[2]
			


# TODO: Propogate Errors to b


	# In[14]:

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
		

	M = sigma_Mach	
	b1 = b_method1(sigma_s,M)
	if(fit_lognormal_M) :
		return sigma_s,error_sigma,M,error_M,b1	
	else :
		return sigma_s,error_sigma,M,b1






if __name__ == "__main__":

	ap = argparse.ArgumentParser(description='Command Line Inputs for b computing script.')
	ap.add_argument('-directory', metavar='Output Directory', action='store', type=str,
		default=os.path.realpath('./')+'/',help='Output files directory')
	ap.add_argument('-N', metavar='Resolution', action='store', type=int,default=200,help='Resolution of Simulation',dest='N')
	ap.add_argument('-output', metavar='Output Directory of table', action='store', type=str,
	                default=os.path.realpath('./')+'/',help='Directory to save table.')
	ap.add_argument('-tstart', metavar='t0', action='store', type=int,default=100,
	            help='Starting timestep.')
	ap.add_argument('-tend', metavar='tend', action='store', type=int,default=300,
	            help='Last timestep.')
	args = vars(ap.parse_args())


	system_t0 = time.time()
	directory = os.path.abspath(args['directory']) + '/'
	output = os.path.abspath(args['output']) + '/'
	N = args['N']
	tstart = args['tstart']
	tfinish = args['tend']
	time = tstart

	SIGMA_S = np.zeros(tfinish-tstart+1)
	ERROR_S = np.zeros(tfinish-tstart+1)
	MACH_NO = np.zeros(tfinish-tstart+1)
	if(fit_lognormal_M) :
		ERROR_M = np.zeros(tfinish-tstart+1)
	B_VALUE = np.zeros(tfinish-tstart+1)
	B_ERROR = np.zeros(tfinish-tstart+1)






	if(use_Density_Threshold) :
	    output = output + 'Density_Threshold/'
	    if not os.path.exists(output):
	        os.makedirs(output)
	else :
	    output = output + 'Region_Threshold/'
	    if not os.path.exists(output):
	        os.makedirs(output)


	while time <=tfinish :

		if(fit_lognormal_M) :
			sigma_s,error_sigma,M,error_M,b1 = calculate_b(directory,time,N)
		else :
			sigma_s,error_sigma,M,b1 = calculate_b(directory,time,N)
			error_Mach = 0.0
		SIGMA_S[time-tstart] = sigma_s
		ERROR_S[time-tstart] = error_sigma
		MACH_NO[time-tstart] = M 
		if(fit_lognormal_M) :
			ERROR_M[time-tstart] = error_Mach
		B_VALUE[time-tstart] = b1 
		if(error_sigma==100.0) :
			B_ERROR = 100.0
		else :	
			B_ERROR[time-tstart] = std_error(sigma_s,M,error_sigma,error_Mach)
		print("%d \n"%time)
		time = time + 1


	## B Values Obtained
	time_column = np.arange((tstart-100)*10,(tfinish-100+1)*10,10)
	if(fit_lognormal_M) :
		header = '1. Time(Kyr) \t sigma_s \t Error_s \t Mach No \t Error_M \t b \t error_b \n\n'
		text_columns = np.column_stack((time_column,SIGMA_S,ERROR_S,MACH_NO,ERROR_M,B_VALUE,B_ERROR))
		np.savetxt(output+'b_timevariation_LognormalM.dat',text_columns,fmt='%0.2f',delimiter='\t',header=header)
	else :
		header = '1. Time(Kyr) \t sigma_s \t Error_s \t Mach No \t b \t error_b \n\n'
		text_columns = np.column_stack((time_column,SIGMA_S,ERROR_S,MACH_NO,B_VALUE,B_ERROR))
		np.savetxt(output+'b_timevariation.dat',text_columns,fmt='%0.2f',delimiter='\t',header=header)	

	print('Time Variation of b Calculated\n\n')	
	import time as time
	system_t1 = time.time()
	total_time = system_t1-system_t0
	print('\n Time Taken = %f seconds'%total_time)



