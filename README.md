## Code to support paper "Leader cells in collective chemotaxis: optimality and tradeoffs"
### Code by Austin Hopkins
These files include the summary data from the paper on Leader Cells in Collective Chemotaxis and code to run the necessary simulations.

AnalyzeData contains a script that will produce the summary data from a given simulation run.

ReplicateFigures contains a jupyter notebook for generating the figures in the paper from the summary data. The summary data that was used in the paper is also provided.

RunSimulation contains a python script that will simulate a linear cluster for a follower noise level delf, a gradient width h, and a size N. It also contains the necessary tables for computing tau and sigma at intermediate leader/follower noise levels.

SingleAngleCorrelationTimes contains a python script that can be used to generate a decade of tau(noise) and sigma(noise). This is how the tables used in the simulation were generated.
