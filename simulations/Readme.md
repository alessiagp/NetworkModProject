# Network Modeling and Simulation final project

This folder contains the scripts used to run the simulation of our ODE system and the scripts used to plot the results. The folder *sim_data* cresults of the simulation results in form of CSV files, whose column represent: time, number of cancer cells, number of effector cells, number of drug molecules. The file names indicate the conditions, i.e.

- No treatment
- Cyt_high: Cyt 62.5 mg/kg 3 days
- Cyt_low: Cyt 0.12 mg/kg 5 days
- Ibr_high: Ibr 9 mg/kg days 1-5 and 8-10
- Ibr_low: Ibr 18 mg/kg days 1-5 and 8-10
- Cyt_Ibr: Cyt 62.5 + Ibr 9 mg/kg days 1-5 and 8-10

These were taken directly from the paper. Additionally we simulated the following conditions

- Ibr_25: Ibr 25 mg/kg days 1-5 and 8-10
- ??

The folder *figures* contains the plots of the CSV files obtained using the script *plotter.py* found in the main folder.

### model_no_trt

This script contains a simplified model without treatment.

### model_trt

This script contains the final model of ODEs, solved using the ODE45 subroutine of Matlab that implements a 4-th order adaptive step Runge-Kutta integration. The options for the numerical solution of the model are: 'MaxStep'=1 to avoid that the integration *skips* a treatment. 'RelTol' = 1e-2 and 'AbsTol' = 1e-4 to avoid numerical errors due to the abrupt change in the number of drug molecules when a treatment is performed.

### sim_data directory

In the files with the *_2m* suffix, from the starting condition $X_0=(5 \times 10^4,2500,0)$, the simulation is run without treatment for 2 weeks of simulated time. Then the treatment is carried out as described above, depending on the drug. In the *figures* folder it is possible to find the plots referring to the different simulations.

The simulations were carried out using both O(1e18) number of molecules or O(1e15) number of molecules, in the latter case the suffix *_1e15drug* was added.

Assuming that the Cytotoxicity rate in the presence of drug $\mu_{AC}$ for Ibr scales as the line fitted to the points sonsidered in the paper, i.e., (9,0.0041) and (18,0.0042), we get: $mu_{AC}=0.00001 \cdot dose + 0.004$.  Then, to obtain $\mu_{AC}=0.01$ we would need a dose of $600 mg/kg$. On the other hand, with a more conservative scaling (i.e. an improvement of 0.001 at each doubing of the concentration) to obtain the same effect we would need a dose of $9^{59}$ mg/kg (which is insane).
in ref [32] and [34] of the main paper, a dose of 25 mg/kg per day of Ibr is used, and results in >90% inhibition of BTK. The $mu_{ac}$ associated with this dose is $\mu_{ac} = 0.0043$ and the number of molecules per dose is $0.5 \times 1.4 \times 10^{18}$.
