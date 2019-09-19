import numpy as np
from netCDF4 import Dataset
import sys

def direct():
    # where the data file was output
    return '../RunSimulation/'
	
def save_direct():
    # where you want to save the summary data
    return './'
	
def get_all_Vcs(Vcses_tot):
    Vc_list = []
    for Vcses in Vcses_tot:
        Vc_list.append(Vcses.compressed())
    return Vc_list
	
def bootstrap_func_mean(X):
	means = []
	for i in range(10000):
		boot = np.random.choice(X,size=len(X))
		means.append(boot.mean())
	means = np.array(means)
	# return 10,000 deltas along with the boostrap mean (for opts purposes)
	return means-X.mean(),means

def get_confidence_mean(X,perc):
	deltas,boot_means = bootstrap_func_mean(X)
	mu = np.mean(X)
	return mu,mu-np.percentile(deltas,100-perc),mu-np.percentile(deltas,perc),boot_means
	
def do_bootstrap_Nl(array,perc,save_directory,filename):
    means = []
    uppers = []
    lowers = []
    for x in array:
        m,l,u,boot_means = get_confidence_mean(x,perc)
        means.append(m)
        uppers.append(u)
        lowers.append(l)
    means = np.array(means)
    uppers = np.array(uppers)
    lowers = np.array(lowers)
    np.save(save_directory+'m_'+filename+'.npy',means)
    np.save(save_directory+'u_'+filename+'.npy',uppers)
    np.save(save_directory+'l_'+filename+'.npy',lowers)
	
def do_bootstrap_Vcses_tot(Vcses_tot,perc_vc,perc_opts,save_directory,filename,type):
	means = []
	uppers = []
	lowers = []
	all_boot_means = []
	for x in Vcses_tot:
		m,l,u,boot_means = get_confidence_mean(x,perc_vc)
		means.append(m)
		uppers.append(u)
		lowers.append(l)
		all_boot_means.append(boot_means)
	means = np.array(means)
	uppers = np.array(uppers)
	lowers = np.array(lowers)
	np.save(save_directory+'m_'+filename+'.npy',means)
	np.save(save_directory+'u_'+filename+'.npy',uppers)
	np.save(save_directory+'l_'+filename+'.npy',lowers)
	
	all_boot_means = np.array(all_boot_means)
	opts = np.where(all_boot_means==np.max(all_boot_means,axis=0))[0]
	np.save(save_directory+'opts'+type+'.npy',opts)

N = int(sys.argv[1])
Nls = np.arange(0,51,1)

directory = direct()
save_directory = save_direct()
		
ds = Dataset(directory + 'data.nc')
		
Vcses = ds['Vx'][:]
Vcses_tot = get_all_Vcs(Vcses)
do_bootstrap_Vcses_tot(Vcses_tot,2.5,16,save_directory,'Vx','Vx')
		
CIs = ds['CI'][:]
CIs_tot = get_all_Vcs(CIs)
do_bootstrap_Vcses_tot(CIs_tot,2.5,16,save_directory,'CI','CI')
		
taus_Nl = ds['taus_1'][:]
do_bootstrap_Nl(taus_Nl,2.5,save_directory,'taus')
do_bootstrap_Nl(1/taus_Nl,2.5,save_directory,'taus_inv')
		
taus2_Nl = ds['taus_2'][:]
do_bootstrap_Nl(taus2_Nl,2.5,save_directory,'taus2')
do_bootstrap_Nl(1/taus2_Nl,2.5,save_directory,'taus2_inv')
		
ds.close()	
