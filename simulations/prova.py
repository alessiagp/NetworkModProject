import numpy as np
import matplotlib.pyplot as plt

#file = "Ibr_Cyt_2m_1e15drug.txt"
file = "Cyt_infusion_low_e0.txt"

data = np.genfromtxt("sim_data/" + file, delimiter=",", skip_header=1)
t = data[::,0]

t = np.asarray(t)

#### First treatment (7*24)
t_1 = t[1344:85850]
diffs = np.asarray(t_1[1:] - t_1[:-1])
print(np.mean(diffs), min(diffs))

#### Second treament (3h infusions)  672 && t < 675

t_2 = t[122109:123653]
diffs_2 = np.asarray(t_2[1:] - t_2[:-1])
print(np.mean(diffs_2), min(diffs_2))