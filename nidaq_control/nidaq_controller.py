if __name__ == '__main__':
    import nidaqmx as nidaq

else:
    from . import nidaqmx as nidaq

import numpy as np


# meine Funktionen
import matplotlib.pyplot as plt

with nidaqmx.Task() as task:
    task.ai_channels.add_ai_voltage_chan("Dev1/ai0")
    task.ai_channels.add_ai_voltage_chan("Dev1/ai1")
    task.timing.cfg_samp_clk_timing(23127, sample_mode=nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan=50000)
    a = np.array(task.read(number_of_samples_per_channel=50000))

ratios = np.zeros(shape=(50000, 2))
for i in range(len(a[0])):
    if a[0,i] > 0.1*a[0].max():
        ratios[i] = np.array([i, a[0,i] / a[1,i]])

plt.plot(a[0], alpha=0.5)
plt.plot(a[1], alpha=.5)
plt.plot(ratios[:,0], ratios[:,1], marker="x", linestyle="")
plt.show()


mean = 0
i = 0
for entry in ratios[:,1]:
    if entry > 0:
        mean += entry
        i += 1
mean = mean / i

std = 0
for entry in ratios[:,1]:
    if entry > 0:
        std += (entry-mean)**2
std = np.sqrt(1/(i-1) * std)

print(mean, std)
