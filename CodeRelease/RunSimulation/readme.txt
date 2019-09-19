The requirements for this script are sys, numpy, scipy and netcdf4-python

netcdf4-python can be installed:
with pi: pip install netcdf4
with conda: conda install -c anaconda netcdf4

run_simulation.py takes the follower noise (in radians), the gradient width parameter h (in cell diameters), and the cluster size (in number of cells) and simulates the cluster in a line geometry as in the paper.

The variable N_configs sets the number of trials per leader number you would like to run. By default, this is set to 5. In the paper, we used 50 to generate statistics for the correlation times.

For example 'python run_simulation.py .314 10 50' will run a simulation of a 50 cell cluster in a line whose followers have a noise of .314 radians in a gradient with h=10. On my machine, this set of parameters takes about 20 seconds per trial. The total time will be (time per trial) x (number of trials per leader number) x (number of leaders).

The files delta_sigma_sq.npy and tau_corr.npy contain data from the simulations of angles relaxing with periodic boundaries. This is used for the numerical interpolation scheme for determining the parameters sigma(uncertainty) tau(uncertainity) so that the cell's polarity correlation time = 1.
