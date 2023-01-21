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

### model_cyt_infusion

This treatment protocols are taken from [https://www.nejm.org/doi/full/10.1056/Nejmoa1010222]


**Cyt_low**: 200 mg * 24h / m^2 for 7 consecutive days + 1000 mg * 3h / m^2 every 12 h for 6 days.

200 mg * 24h / m^2 = 66.7 mg * 24h / kg = 2.8 mg * h / kg, this corresponds to 1.4 mg/mouse in a 24h infusion, which corresponss to 1.4 * 2.4e18 molecules/mouse per one dose, since the treatments lasts one week, the total numbe of molecules will be: 9.8 * 2.4e18 : $mu_{AC} = 0.00148$

1000 mg/m^2 = 333 mg * 3h / kg = 111 mg * h / kg, this corresponds to 2.3 mg/mouse in a 3h infusion, which corresponds to 2.3 * 2.4e18 molecules/mouse per dose: $mu_{AC} = 0.021$

**Cyt_high**: 1000 mg * 3h / m^2 every 12 h for 5 days + 2000 mg * 6h / m^2 every 12 h on days 12, 13, 15, 17.

1000 mg/m^2 = 333 mg * 3h / kg = 111 mg * h / kg, this corresponds to 2.3 mg/mouse in a 3h infusion, which corresponds to 2.3 * 2.4e18 molecules/mouse per dose: $mu_{AC} = 0.021$

2000 mg/m^2 = 667 mg * 6h / kg = 111 mg * / kg, this corresponds to 4.6 mg/mouse in a 6h infusion, which corresponds to 4.6 * 2.4e18 molecules/mouse per dose: $mu_{AC} = 0.021$

To estimate $mu_{AC}$ in these conditions, we used the line $mu_{AC} = 0.00018 * dose [mg * h/kg] + 0.00098$.

To calculate the number of molecules per time-step dt, we used the mean value of the time-step when the treatment asre delivered. This was calculated by looking at the vecors of time of some "dummy" simulations.

### sim_data directory

In the files with the *_2m* suffix, from the starting condition $X_0=(5 \times 10^4,2500,0)$, the simulation is run without treatment for 2 weeks of simulated time. Then the treatment is carried out as described above, depending on the drug. In the *figures* folder it is possible to find the plots referring to the different simulations.

The simulations were carried out using both O(1e18) number of molecules or O(1e15) number of molecules, in the latter case the suffix *_1e15drug* was added.

Assuming that the Cytotoxicity rate in the presence of drug $\mu_{AC}$ for Ibr scales as the line fitted to the points sonsidered in the paper, i.e., (9,0.0041) and (18,0.0042), we get: $mu_{AC}=0.00001 \cdot dose + 0.004$.  Then, to obtain $\mu_{AC}=0.01$ we would need a dose of $600 mg/kg$. On the other hand, with a more conservative scaling (i.e. an improvement of 0.001 at each doubing of the concentration) to obtain the same effect we would need a dose of $9^{59}$ mg/kg (which is insane).
in ref [32] and [34] of the main paper, a dose of 25 mg/kg per day of Ibr is used, and results in >90% inhibition of BTK. The $mu_{ac}$ associated with this dose is $\mu_{ac} = 0.0043$ and the number of molecules per dose is $0.5 \times 1.4 \times 10^{18}$.

In the files with the *_e0* suffix, a term for the natural production of effector cells was addes. This term $c0 = 4.2 \ [cells/h]$ was derived from ref. [36].


Cyt_low = 200 mg / m^2 = 66.7 mg * 24h / kg = 2.8 mg * h / kg	Cyt_low = 1.4 mg / mouse (24h infusion)
Cyt_high = 2500 mg/m^2 = 833 mg * 3h / kg = 277 mg * h / kg 3h - 2times a day Cyt_high = 17.5 mg / mouse


