# Simulations

This folder contains the results of the simulations in form of CSV files. The columns are: time, number of cancer cells, number of effector cells, number of drug molecules. The file names indicate the conditions, i.e.

- Cyt_high: Cyt 62.5 mg/kg 3 days
- Cyt_low: Cyt 0.12 mg/kg 5 days
- Ibr_high: Ibr 9 mg/kg days 1-5 and 8-10
- Ibr_low: Ibr 18 mg/kg days 1-5 and 8-10
- Cyt_Ibr: Cyt 62.5 + Ibr 9 mg/kg 8 days

These were taken from the paper.
In the files with the *_2months* suffix, from the starting condition of X0=(5e4,2500,0), the simulation is run without treatment for 3 weeks of simulated time. Then the treatment is carried out as described before.
In the *figures* folder it is possible to find the plots referring to the different simulations. The starting condition is X0=(5e4,2500,0) in all the cases.
