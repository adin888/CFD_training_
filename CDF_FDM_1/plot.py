import matplotlib.pyplot as plt
import numpy as np

[q1,q2,q3,q4] = np.loadtxt("C:/Users/Administrator/source/repos/CFD/CDF_FDM_1/field_final.csv",delimiter=",",skiprows=1,unpack=True)

plt.plot(q1, q2)
plt.show()
