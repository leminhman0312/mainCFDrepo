import matplotlib.pyplot as plt
import numpy as np

# set range for Pe
Pe_vals = [-50,-5,1e-6,5,50]

# assign x values
x = np.linspace(0,1)

# set TP and TE
TP = 1
TE = 0

for Pe in Pe_vals:
    T = (TE-TP)*(np.exp(Pe*(x-x[0])/(x[-1]-x[0]))-1) \
        /(np.exp(Pe)-1)
    plt.plot(x, T, label=str(int(Pe)))

plt.xlabel(r"$x/\Delta x$")
plt.ylabel(r"$T$")
plt.legend()
plt.savefig("pic/pectet_visual_pipe.png")
