### Extra-Pair Mating and the Evolution of Cooperative Behaviours
### Agnieszka Rumińska, Christian Jørgensen, Sigrunn Eliassen
### 15.11.2022
### file 3/3

import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D      

"""Quickly depickling *details.pickle files.

Structure of the pickle file:
carea : list
	all c_1 values in the population, extracted from males
careb : list
	all c_2 values in the population, extracted from males
carea_f : list
	all c_1 values in the population, extracted from females
careb_f : list
	all c_2 values in the population, extracted from females
defa : list
	all d_1 values in the population, extracted from males
defb : list
	all d_2 values in the population, extracted from males
defc : list
	all d_3 values in the population, extracted from males
defa_f : list
	all d_1 values in the population, extracted from females
defb_f : list
	all d_2 values in the population, extracted from females
defc_f : list
	all d_3 values in the population, extracted from females
epc_f : list
	all epc values in the population, extracted from females
epc_m : list
	all epc values in the population, extracted from males
		"""

def check_traits():
	plt.figure(figsize=(6, 10))
	plt.subplots_adjust(left=0.135, bottom=0.1, right=None, top=0.9, hspace=0.4)
	plt.suptitle('Male care and territorial defence', size=18)

	plt.subplot(311)
	plt.plot(time, c_1_means, marker = 'o', ms=1, color='C0')
	plt.plot(time, c_2_means, marker = 'o', ms=1, color='C1')
	legend = [Line2D([0], [0], color='C0', lw=4), Line2D([0], [0], color='C1', lw=4)]   
	plt.legend(legend, ['$c_{1}$', '$c_{2}$'], prop={'size': 12}, loc='lower left', ncol=1)
	plt.grid(axis='y')
	plt.ylabel('genetic value', size=16)
	plt.xlabel('time ($10^{2}$ years)', size=16)

	plt.subplot(312)
	plt.plot(time, d_1_means, marker = 'o', ms=1, color='C2')
	plt.plot(time, d_2_means, marker = 'o', ms=1, color='C3')
	legend = [Line2D([0], [0], color='C2', lw=4), Line2D([0], [0], color='C3', lw=4)]   
	plt.legend(legend, ['$d_{1}$', '$d_{2}$'], prop={'size': 12}, loc='lower left', ncol=1)
	plt.grid(axis='y')
	plt.ylabel('genetic value', size=16)
	plt.xlabel('time ($10^{2}$ years)', size=16)

	plt.subplot(313)
	plt.plot(time, d_3_means, marker = 'o', ms=1, color='C4')
	legend = [Line2D([0], [0], color='C4', lw=4)]   
	plt.legend(legend, ['$d_{3}$'], prop={'size': 12}, loc='lower left', ncol=1)
	plt.grid(axis='y')
	plt.ylabel('genetic value', size=16)
	plt.xlabel('time ($10^{2}$ years)', size=16)

	plt.savefig('EPP_genes.png', dpi=300)
	plt.show()

breeding_pairs_number = 24*220
carea, careb, defa, defb, defc = [], [], [], [], []
const = 0

for i in range(0, len(sys.argv)-1):
	file = open(sys.argv[i+1], "rb")
	print(sys.argv[i+1][:])
	while True:
		try:
			statistics = (pickle.load(file))
			carea.append(statistics[0])
			careb.append(statistics[1])
			defa.append(statistics[4])
			defb.append(statistics[5])
			defc.append(statistics[6])
			const+=1
		except EOFError:
			break
	file.close()


c_1_means, c_2_means = [], []
d_1_means, d_2_means, d_3_means = [], [], []
for a in range(0, len(carea)):
	c_1_means.append(np.mean(carea[a]))
	c_2_means.append(np.mean(careb[a]))
	d_1_means.append(np.mean(defa[a]))
	d_2_means.append(np.mean(defb[a]))
	d_3_means.append(np.mean(defc[a]))

if len(c_1_means) == const:
	time = np.arange(0, const, 1)
	check_traits()
else:
	print('Error')