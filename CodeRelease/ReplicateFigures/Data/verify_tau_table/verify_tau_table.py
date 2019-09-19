import numpy as np
import time
import sys
from scipy.optimize import curve_fit

def o_u(sigma,tau,dt):
    tmax = int(3000/dt)
    thetas = []
    theta = np.random.uniform(-np.pi,np.pi)
    thetas.append(theta)
    ts = [0]
    for i in range(tmax):
        theta += -((((theta+np.pi)%(2*np.pi))-np.pi)/tau)*dt + sigma*np.sqrt(dt)*np.random.normal(0,1)
        thetas.append(theta)
        ts.append((i+1)*dt)
    return np.array(ts),np.array(thetas)

def vel_correlation(vcor_samples,Vcs):
    vcor = []
    msq = np.dot(np.mean(Vcs,axis=0),np.mean(Vcs,axis=0))
    vcor.append(np.mean(Vcs[:,0]*Vcs[:,0]+Vcs[:,1]*Vcs[:,1])-msq)
    for i in vcor_samples[1:]:
        vcor.append(np.mean((Vcs[:-i:,0]*Vcs[i::,0])+(Vcs[:-i:,1]*Vcs[i::,1]))-msq)
    return np.array(vcor)
	
def tau_corr_sim(delta_squared,trials):
    dt = 0.01
    Nsamples = 600
    sigma = sig_func(delta_squared)
    tau = tau_func(delta_squared)
    vcor_samples = np.append(np.array([0]),np.unique(np.logspace(0,np.log10(Nsamples),num=int(Nsamples)/10,dtype=int)))
    samps = []
    for i in range(trials):
        ts,thetas = o_u(sigma,tau,dt)
        Vcs = np.array(np.array([np.cos(thetas),np.sin(thetas)]).T)
        vcorr = vel_correlation(vcor_samples,Vcs[int(3*tau/dt):,:])
        p,cov = curve_fit(exp_func,ts[vcor_samples],vcorr/vcorr[0])
        samps.append(p[0])
    samps = np.array(samps)
    return samps.mean(),samps.std()/np.sqrt(trials)
	
def exp_func(t,tau):
    return np.exp(-t/tau)
	
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

known_sigmas = np.load('/scratch/groups/bcamley1/data/tau_table/delta_sigma_sq.npy')
known_taus = 1/np.load('/scratch/groups/bcamley1/data/tau_table/tau_corr.npy')
	
deltas_squared = np.logspace(-3,2,num=100)

tcorrs = []
tcorrs_err = []
for delta_sq in deltas_squared:
    s = time.time()
    tcorr, err = tau_corr_sim(delta_sq,100)
    tcorrs.append(tcorr)
    tcorrs_err.append(err)
    print(delta_sq,tcorr,err,s-time.time())
	
np.save('./deltas_squared_verify.npy',deltas_squared)
np.save('./tcorrs_verify.npy',tcorrs)
np.save('./tcorrs_err_verify.npy',tcorrs_err)
