# top-geom-meas
Code and data used in the paper "Topological transition in measurement-induced geometric phases"

The description of files
-------------------------

1. geophase_mc.C - a program in C for performing Monte Carlo simulations and saving the dependence of the averaged phase on the polar angle (theta)

2. protocol1_cXXX.dat - the dependence of <exp(i \chi_{geom})> on theta at the measurement strength parameter c=XXX. The first column is theta, the second is arg(<exp(i \chi_{geom})>)/pi, the third is abs(<exp(i \chi_{geom})>)

3. protocol2_cXXX.dat - the dependence of <exp(2 i \chi_{geom})> on theta at the measurement strength parameter c=XXX. The first column is theta, the second is 0.5\*arg(<exp(2 i \chi_{geom})>)/pi, the third is abs(<exp(2 i \chi_{geom})>)

4. compute_phases_histogram.py - a program in python that performs the Monte Carlo simulation for a fixed polar angle (theta) and different measurement strengths, and outputs the averaged geometric phases and the visibility, as well as the histogram of the phase distribution, needed to produce Fig. 3. 

5. plotting_phasespaper_pi_4.txt - Matrix with information of the simulated average geometric phase. Each row corresponds to the simulation of a different measurement strength. Each row has the following structure: 
(measurement strength, <\cos(2 \chi_{geom})>, std of \cos(2 \chi_{geom}), <\sin(2 \chi_{geom})>, std of \sin(2 \chi_{geom}), arg(<exp(2 i \chi_{geom})>)/2, average probability of succesfull final projection, average prob of succ of final projection  times abs( <exp(2 i \chi_{geom})>) ). < > are averages wheighted with the probability of a successful final projection. 

6. whistogram_paperpi_4 - data for the histogram of Fig. 3: Each row corresponds to a different measurement strength (see plottling_phasespaper_pi_4.txt). Each row has 180 entries (resolution), encoding the number of phases simulated in the corresponding interval, and wheighted by the probability of a succesful final measurement.

7. plotfig3.nb - mathematica file to produce Fig. 3: Plot of the geometric phase of the post-selected tajectory and the average geometric phase together with the histogram of the phase distribution. Produces also the inserts: Legend of the histogram and plot of the probability of the post-selected phase and the visibility of the phase. theta = Pi/4


Authorship statement
-------------------------
Files 1, 2, 3 have been produced by Thomas Wellens. Files 4, 5, 6, 7 have been produced by Valentin Gebahrt. 
