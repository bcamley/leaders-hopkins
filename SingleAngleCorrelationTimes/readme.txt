generate_tau_table.py computes the polarity correlation time Tp as a function of angular uncertainity at a given parameter tau. From Tp(uncertainity,tau), the correct parameters sigma(uncertainity) and tau(uncertainity) can be determined so that Tp=1 at all levels of uncertainity.

generate_tau_table.py expects n as an argument, which specifies the order of magnitude over which you wish to compute Tp. That is 'python generate_tau_table.py -1' will compute Tp(uncertainity,tau) for uncertanties between 10^-1 and 10^0. generate_tau_table.py assumes the parameter tau=1 because Tp/tau has a universal scaling relation with uncertainty = sigma^2*tau/2. For a given n, the grid spacing will be 10^(n-1).

This script was used to generate the values in tau_corr.npy that were used in the interpolation scheme for the cluster simulations.
