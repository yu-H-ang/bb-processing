import numpy as np
import matplotlib.pyplot as plt
import pywt

with open("output.txt", "r") as f:
    tmp = [float(i) for i in f.read().split()]
kn = 160
ts = 338
te = 3972
time_length = te - ts + 1
data = np.zeros((time_length,kn))

for i in range(time_length):
    for j in range(kn):
        data[i,j] = tmp[i*kn+j]
data = np.transpose(data)

plt.imshow(data[:,500-ts:], cmap='Reds', origin='lower', aspect='auto', extent=(500,te,-12.6,12.6), clim=[0,1])
plt.colorbar()
plt.title("averaged volume fraction")
plt.xlabel("time(ms)")
plt.ylabel("vertical coordinate(mm)")
plt.show()
