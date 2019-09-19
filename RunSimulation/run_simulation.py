import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
import sys
from netCDF4 import Dataset
import time
	
def C(x,h,Cmax,Cmin):
    return (Cmax-Cmin)*(1/2)*(1+np.tanh(x/h))+Cmin

def Cphi_g(phi_r,x,h,Cmax,Cmin):
    return (Cmax-Cmin)*(1/2)*(1+np.tanh((x+(1/2)*np.cos(phi_r))/h))+Cmin

def DcDphi_g(phi_r,x,h,Cmax,Cmin):
    return -(Cmax-Cmin)*(np.cosh((x+(1/2)*np.cos(phi_r))/h)**(-2))*np.sin(phi_r)/(4*h)

def FI_phi_r(phi_r,x,h,Cmax,Cmin,kd):
    con = Cphi_g(phi_r,x,h,Cmax,Cmin)
    dcon = DcDphi_g(phi_r,x,h,Cmax,Cmin)
    return (kd*dcon**2)/(con*(con+kd)**2)

def sig_deltas(x,h,Cmax,Cmin):
	FIs = []
	for i in x:
		FIs.append(Nr*integrate.quad(FI_phi_r, 0, 2*np.pi, args=(i,h,Cmax,Cmin,kd))[0]/(2*np.pi))
	FIs = np.array(FIs,dtype=float)
	return 1/np.sqrt(FIs)

def dtheta_f(Pfs,Vc,tau):
    dot = np.einsum('i,ji->j',Vc,Pfs)
    det = Vc[0]*Pfs[:,1]-Vc[1]*Pfs[:,0]
    return -np.arctan2(det,dot)/tau
	
def vel_correlation(Vcs):
	vcor = []
	msq = np.dot(np.mean(Vcs,axis=0),np.mean(Vcs,axis=0))
	vcor.append(np.mean(Vcs[:,0]*Vcs[:,0]+Vcs[:,1]*Vcs[:,1])-msq)
	for i in vcor_samples[1:]:
		vcor.append(np.mean((Vcs[:-i:,0]*Vcs[i::,0])+(Vcs[:-i:,1]*Vcs[i::,1]))-msq)
	return np.array(vcor)

def exp_func(x,tau):
    return np.exp(-x/tau)
	
def get_tau2(vcor):
	# vcor value and time just above 1/2
	v1 = vcor[np.where(vcor/vcor[0]>=1/2)[0][-1]]/vcor[0]
	t1 = vcor_samples[np.where(vcor/vcor[0]>=1/2)[0][-1]]*delta_t
	# vcor value and time just below 1/2
	v2 = vcor[np.where(vcor/vcor[0]<1/2)[0][0]]/vcor[0]
	t2 = vcor_samples[np.where(vcor/vcor[0]<1/2)[0][0]]*delta_t
	# slope of line between points
	slope = (v2-v1)/(t2-t1)
	# find where interpolated line equals 1/2
	tfin = ((1/2)-v1)/slope+t1
	# return rescaled value to 1/e
	return tfin*.5*np.exp(1)
	
def tau_func(del_sq):
	# careful - input needs to be a float
    conds = [del_sq < 0.01, np.logical_and(del_sq >= 0.01,del_sq <= 10), del_sq > 10]
    funcs = [lambda x: 1,lambda x: np.interp(x,known_sigmas,known_taus),lambda x: x]
    return np.piecewise(del_sq,conds,funcs)
	
def sig_func(del_sq):
	# careful - input needs to be a float
    conds = [del_sq < 0.01, np.logical_and(del_sq >= 0.01,del_sq <= 10), del_sq > 10]
    funcs = [lambda x: np.sqrt(2*x),lambda x: np.sqrt(2*x/np.interp(x,known_sigmas,known_taus)),lambda x: np.sqrt(2)]
    return np.piecewise(del_sq,conds,funcs)

def run_simulation(Nl,N,N_step,Cmax,Cmin,h,dt,delnu):
	# Initialize
	tau = 1

	total_steps = round(N_step*tau/dt)

	potential_leaders = np.linspace(0,-N+1,num=N)
	followers = np.linspace(Nl,N-1,num=N-Nl)

	# initialize leader sigmas and taus
	potential_ldeltas = sig_deltas(potential_leaders,h,Cmax,Cmin)
	ldeltas = potential_ldeltas[np.argsort(potential_ldeltas)][:Nl]
	
	lsigs = sig_func(ldeltas**2)
	taus_l = tau_func(ldeltas**2)

	# initialize follower sigma and tau
	fsig = sig_func(delnu**2)
	tau_f = tau_func(delnu**2)


	# initialize leader velocities
	thetas_l = np.random.uniform(-np.pi,np.pi,Nl)
	Pls = np.array([np.cos(thetas_l),np.sin(thetas_l)]).T

	# initialize follower velocities
	thetas_f = np.random.uniform(-np.pi,np.pi,N-Nl)
	Pfs = np.array([np.cos(thetas_f),np.sin(thetas_f)]).T

	# initial cluster velocity
	Vc = (1/N)*(np.sum(Pls,axis=0)+np.sum(Pfs,axis=0))

	Vcs = []
	Vcs.append(Vc)

	# generate all noises for run
	leader_noise = lsigs[np.newaxis,:]*np.sqrt(dt)*np.random.normal(0,1,(total_steps-1,Nl))
	follower_noise = fsig*np.sqrt(dt)*np.random.normal(0,1,(total_steps-1,N-Nl))

	# run a time step
	i = 1
	while i < total_steps:
		# good leaders update angle
		thetas_l += -(thetas_l/taus_l)*dt + leader_noise[i-1]
		thetas_l = ((thetas_l+np.pi)%(2*np.pi)-np.pi)
		Pls = np.array([np.cos(thetas_l),np.sin(thetas_l)]).T

		# followers update their angle
		thetas_f += dtheta_f(Pfs,Vc,tau_f)*dt + follower_noise[i-1]
		Pfs = np.array([np.cos(thetas_f),np.sin(thetas_f)]).T
				
		Vc = (1/N)*(np.sum(Pls,axis=0)+np.sum(Pfs,axis=0))
		
		Vcs.append(Vc)
		
		i += 1

	Vcs = np.array(Vcs)
	return np.linspace(0,N_step,num=total_steps),Vcs

# Input and labeling
delnu = float(sys.argv[1])
h = float(sys.argv[2])
N = int(sys.argv[3])
slabel = str(int(100*delnu/(2*np.pi)))
# number of times you want to run the simulation per leader number
N_configs = 5

#### Description of output
print('Cluser Size N =',N)
print('follower noise delnu =',delnu)
print('gradient width h =',h)
print("""--------------------------------------------------------------------
********************************************************************
--------------------------------------------------------------------
""")
Vcses_Nl = []
CIs_Nl = []
Nls = np.arange(0,N+1,1)
Nsamples = 60000
vcor_samples = np.append(np.array([0]),np.unique(np.logspace(0,np.log10(Nsamples),num=int(Nsamples/100),dtype=int)))
delta_t = 0.01
Nr = 70000
kd = 1

# Create NETCDF File
fname = './data.nc'
ds = Dataset(fname,'w',data_model='NETCDF4',file_format='NETCDF4')
ds.createDimension("leader",N+1)
ds.createDimension("sample",None)
ds.createDimension("tau_samples",N_configs)
Vx = ds.createVariable("Vx",'f4',("leader","sample"),fill_value=-2)
CI = ds.createVariable("CI",'f4',("leader","sample"),fill_value=-2)
T1 = ds.createVariable("taus_1",'f4',("leader","tau_samples"))
T2 = ds.createVariable("taus_2",'f4',("leader","tau_samples"))

known_sigmas = np.load('./delta_sigma_sq.npy')
known_taus = 1/np.load('./tau_corr.npy')

for Nl in Nls:
        print('')
        print('Number of leaders is ' + str(Nl))
        print('-------------------------------------------')
        print('Trial number \t tau (from exponential fit) \t tau (from rescaling)')
        Vcses = []
        CIs = []
        s = 0
        for config in range(N_configs):
                start = time.time()
                ts,Vcs = run_simulation(Nl,N,5200,2,1e-14,h,0.01,delnu)
                vcor = vel_correlation(Vcs[20000:])
                tau = curve_fit(exp_func,vcor_samples*delta_t,vcor/vcor[0])[0][0]
                tau2 = get_tau2(vcor)
                Vc_steady = Vcs[20000:]
                tot_time_steps = len(Vc_steady)
                tot_time = tot_time_steps*delta_t
                intervals = int(np.floor(tot_time/(3*tau)))
                tau_steps = int(np.floor(3*tau/delta_t))
                tmp_vc_list = []
                tmp_ci_list = []
                for i in range(intervals):
                    Vcss = Vc_steady[i*tau_steps:(i+1)*tau_steps,0].mean(axis=0)
                    CIss = Vc_steady[i*tau_steps:(i+1)*tau_steps,0].mean()/(np.sqrt(Vc_steady[i*tau_steps:(i+1)*tau_steps,0]**2+Vc_steady[i*tau_steps:(i+1)*tau_steps,1]**2).mean())
                    tmp_vc_list.append(Vcss)
                    tmp_ci_list.append(CIss)
                Vx[Nl,s:s+intervals] = np.array(tmp_vc_list)
                CI[Nl,s:s+intervals] = np.array(tmp_ci_list)
                T1[Nl,config] = tau
                T2[Nl,config] = tau2
                s += intervals
                print('trial',config,'\t',tau,'\t',tau2)
                print(time.time()-start)
