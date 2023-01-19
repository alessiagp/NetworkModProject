import numpy as np
import matplotlib.pyplot as plt

#file = "Ibr_Cyt_2m_1e15drug.txt"
file = "Ibr_25mg_2m_1e15drug.txt"

data = np.genfromtxt("sim_data/" + file, delimiter=",", skip_header=1)
t = data[::,0]
A = np.around(data[::,1],0)
E = np.around(data[::,2],0)
C = np.around(data[::,3] / 1e8,0)

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(111)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)    
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('Count')
name = file[:-4]
ax1.plot(t,A,lw=1,alpha=0.8,color="darkmagenta",label="A20 cells")
ax1.plot(t,E,lw=1,alpha=0.8,color="royalblue",label="Effector cells")
ax1.plot(t,C,lw=1,alpha=0.8,color="chocolate",label="Drug molecules / 1e8")
ax1.legend(loc='best')
plt.savefig("figures/" + name + ".png")
plt.show()