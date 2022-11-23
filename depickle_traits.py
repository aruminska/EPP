### Extra-Pair Mating and the Evolution of Cooperative Behaviours
### Agnieszka Rumińska, Christian Jørgensen, Sigrunn Eliassen
### 15.11.2022
### file 2/3

import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D     

"""Quickly depickling *traits.pickle files.

Structure of the pickle file:
mean_epc : list
	contains a single value (EPC mean) for each year
sd_epc : list
	contains a single value (standard deviation of EPC) for each year
mean_c : list
	contains a single value (mean expressed male care) for each year
sd_c : list
	contains a single value (standard deviation of male care) for each year
mean_d : list
	contains a single value (mean expressed male defence) for each year
sd_d : list
	contains a single value (standard deviation of male defence) for each year
"""

def check_traits():
	plt.figure(figsize=(10, 6))
	plt.subplots_adjust(left=None, bottom=0.19, right=None, top=0.9, hspace=0.4)

	legend1 = [Line2D([0], [0], color='FireBrick', lw=4),
			Line2D([0], [0], color='LimeGreen', lw=4),
			Line2D([0], [0], color='SteelBlue', lw=4)]  

	plt.suptitle('Traits and genes in the population', size=18)
	plt.plot(time, mean_d, marker = 'o', ms=1, color='SteelBlue')
	plt.plot(time, mean_c, marker = 'o', ms=1, color='LimeGreen')
	plt.plot(time, mean_epc, marker = 'o', ms=1, color='FireBrick')
	plt.legend(legend1, ['epc', 'care', 'defence'], prop={'size': 12}, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3)
	plt.grid(axis='y')
	plt.ylabel('trait value', size=16)
	plt.xlabel('time (years)', size=16)
	plt.savefig('EPP_traits.png', dpi=300)
	plt.show()

breeding_pairs_number = 24*220
mean_epc, sd_epc, mean_c, sd_c, mean_d, sd_d = [], [], [], [], [], []
const = 0

for i in range(0, len(sys.argv)-1):
	file = open(sys.argv[i+1], "rb")
	print(sys.argv[i+1][:])
	while True:
		try:
			statistics = (pickle.load(file))
			mean_epc.append(statistics[0])
			sd_epc.append(statistics[1])
			mean_c.append(statistics[2])
			sd_c.append(statistics[3])
			mean_d.append(statistics[4])
			sd_d.append(statistics[5])
			const+=1
		except EOFError:
			break
	file.close()
time = np.arange(0, const, 1)

check_traits()
