### Extra-Pair Mating and the Evolution of Cooperative Behaviours
### Agnieszka Rumińska, Christian Jørgensen, Sigrunn Eliassen
### 13.10.2022
### file 1/3

# LIBRARIES
import matplotlib.pyplot as plt 
import math
import networkx as nx
import numpy as np 
import numpy.random as nrand
import pickle
import portion 
import statistics as sta
import sys
import time

# CLASSES
class female():
	"""A class that represents female individuals 
	
	Attributes:
	nest_number: int
		position in a social network (social network == neighbourhood == tree)
		females and males occupying the same nest are referred to as a social pair
	tree_number: int
		identifies the social network number within the population
	global_nest_number: int
		global position in the population that consists of many social networks
	neighbours: list
		a list of neighbours that have connections with the individual
		indicates potential extra-pair mates
	eggs: float
		egg production cost; the first component of a female's reproductive investment
	care: float
		care for the fledglings; the second component of a female's reproductive investment
	epc: float
		extra-pair copulation rate (EPC) of a female
		evolving gene, expressed only in females
		denoted in the manuscript by 'x'
	fitness: float
		fitness value returned by W_function()
		used for ranking females during a breeding period
	m_care_a: 
		slope of the male care reaction norm (c_1)
		evolving gene, not expressed in females, but present in their genome
	m_care_b: 
		intercept of the male care reaction norm (c_2)
		evolving gene, not expressed in females, but present in their genome
	m_defence_a: 
		slope in the X direction of the male defence reaction norm (d_1)
		evolving gene, not expressed in females, but present in their genome
	m_defence_b: 
		slope in the Y direction of the male defence reaction norm (d_2)
		evolving gene, not expressed in females, but present in their genome
	m_defence_c: 
		third coefficient of the male defence reaction norm (d_3)
		evolving gene, not expressed in females, but present in their genome
	territory: float
		size of a territory secured by a social partner
	mother: int
		global nest number of an individual's mother 
	father: int
		global nest number of an individual's father
	offspring_list: list of lists
		each sub-list contains female fitness and EPC values, breeding season number and global_nest_number of the father
		the above values are appended separately for each offspring
	breeding_seasons: int
		age of the female, as a number of breeding seasons
	"""
	
	def __init__(self):
		self.care = female_care
		self.epc = starting_epc 
		self.m_care_a, self.m_care_b, = 0., 0.5
		self.m_defence_a, self.m_defence_b, self.m_defence_c = 0, 0, 0.5
		self.mother = -1
		self.father = -1
		self.offspring_list = []
		self.eggs = eggs
		self.breeding_seasons = 1
		self.territory = 0

	def set_global_number(self, num): 
		self.global_nest_number = num
	def get_global_number(self):
		return self.global_nest_number
	def set_tree(self, num):
		self.tree_number = num
	def set_nest_number(self, num):
		self.nest_num = num
	def get_nest_number(self):
		return self.nest_num
	def set_neighbours(self, nlist):
		self.neighbours = []
		for n in nlist:
			self.neighbours.append(self.tree_number*individuals+n)
	def get_neighbours(self):
		return self.neighbours
	def parents(self, mother, father):
		"""Passes parameters to other functions that facilitate genetic inheritance.

			Parameters:
				mother : num
					global_nest_number of the mother
				father : num
					global_nest_number of the mother
		"""
		
		self.mother, self.father = mother, father
		self.set_epc(female_parents[mother].get_epc(), male_parents[father].get_female_epc())
		self.set_male_defence(female_parents[mother].get_male_defence(), male_parents[father].get_defence_params())
		self.set_male_care(female_parents[mother].get_male_care(), male_parents[father].get_care_params())
	def set_epc(self, val1, val2):
		"""Calculates genetic value for x (EPC) from a normal distribution with variance = 0.05.

		Parameters:
			val1 : float
				genetic value of x inherited from mother
			val2 : float
				genetic value of x inherited from father
		Returns: 
			epc : float
		"""

		mean = (val1+val2)/2.
		self.epc = nrand.normal(mean, 0.05)
		if self.epc < 0.:
			self.epc = 0.
		if self.epc > 1.:
			self.epc = 1.
	def get_epc(self):
		return self.epc
	def get_expressed_c(self):
		return self.care
	def get_territory(self):
		if self.territory == 0:
			self.territory = male_parents[self.global_nest_number].calculate_territory()
		return self.territory
	def W_function(self):
		"""Calculates female fitness.
		The function is used for assessing reproductive investment of both parents.
		Finds chances of producing and raising offspring that survive to the next breeding season.

		Returns: 
			fitness : float
		"""

		self.fitness = (self.get_territory()**alpha) \
			* ( self.S_inseason() * self.get_expressed_c()**gamma \
			+ male_parents[self.global_nest_number].S_inseason() * male_parents[self.global_nest_number].get_expressed_c()**gamma )
		return self.fitness
	def get_fitness(self, a, b):
		"""Returns normalised fitness value (between 0.0 and 1.0, used for optimising simulation speed)
		
		Parameters:
			a : float
				Smallest fitness value expressed in the whole population 
			b : float
				Greatest fitness value expressed in the whole population 
		Returns: 
			normalised fitness
		"""

		if a != b:
			return (self.fitness - a)/(b-a)
		else:
			return self.fitness
	def set_male_defence(self, list1, list2): 
		"""Calculates genetic coefficients for the reaction norm (d_1, d_2, d_3).
		New values are drawn from a normal distribution with variance = 0.05.

		Parameters:
			list1 : list
				genetic values of d_1, d_2 and d_3 inherited from mother
			list2 : list
				genetic values of d_1, d_2 and d_3 inherited from father
		Returns: 
			m_defence_a : float
				new slope in the X direction of the male defence reaction norm (d_1)
			m_defence_b : float
				new slope in the Y direction of the male defence reaction norm (d_2)
			m_defence_c : float
				new intercept of the male defence reaction norm (d_3)
		"""

		par1, par2 = population_epc_mean, (1.-population_epc_mean)/nmax
		ang1 = math.atan(list1[0])
		ang2 = math.atan(list2[0])
		mean_a = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_a, 0.05)
		self.m_defence_a = math.tan(new_angle)
		ang1 = math.atan(list1[1])
		ang2 = math.atan(list2[1])
		mean_b = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_b, 0.05)
		self.m_defence_b = math.tan(new_angle)
		Z_val1 = list1[0]*par1 + list1[1]*par2 + list1[2]
		Z_val2 = list2[0]*par1 + list2[1]*par2 + list2[2]
		Z_val_offspring = (Z_val1+Z_val2)/2.
		new_c = Z_val_offspring - self.m_defence_a*par1 - self.m_defence_b*par2
		self.m_defence_c = np.random.normal(new_c, 0.05)
	def get_male_defence(self): 
		return [self.m_defence_a, self.m_defence_b, self.m_defence_c]
	def set_male_care(self, list1, list2):
		"""Calculates genetic coefficients for the reaction norm (c_1, c_2,).
		New values are drawn from a normal distribution with variance = 0.05.

		Parameters:
			list1 : list
				genetic values of c_1 and c_2 inherited from mother
			list2 : list
				genetic values of c_1 and c_2 inherited from father
		Returns: 
			m_care_a : float
				new slope of the male care reaction norm (c_1)
			m_care_b : float
				new intercept of the male care reaction norm (c_2)				
		"""

		par1 = population_epc_mean
		ang1 = math.atan(list1[0])
		ang2 = math.atan(list2[0])
		mean_a = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_a, 0.05)
		self.m_care_a = math.tan(new_angle)
		Y_parent1 = list1[0]*par1 + list1[1]
		Y_parent2 = list2[0]*par1 + list2[1]
		Y_offspring = (Y_parent1+Y_parent2)/2.
		new_b = Y_offspring - self.m_care_a*par1 
		self.m_care_b = np.random.normal(new_b, 0.05)
	def get_male_care(self):
		return [self.m_care_a, self.m_care_b]
	def choose_mate(self):
		"""Calculates the chances of siring an offspring among all males in the neighbourhood.
		Male from the social pair has 1-x chances of siring offspring. The rest (x) is equally distributed among the neighbours.
		The function draws a random number and finds a sire.

		Returns:
			mate : int
				global_nest_number of the chosen mate
		"""

		self.candidates = [[self.global_nest_number, 1 - self.get_epc()]]
		for i in self.neighbours:
			self.candidates.append([i, self.candidates[-1][1] + self.get_epc()/len(self.neighbours)])
		rand = nrand.uniform(0., self.candidates[-1][1])
		for j in range(0, len(self.candidates)):
			if rand <= self.candidates[j][1]:
				mate = self.candidates[j][0]
				break
		return mate 
	def S_function(self):
		"""Calculates between-seasons survival probability of a female.
		If a female invests too much in care and egg production, it results in elevated mortality risk due to energy loss.
		S_function is used to assess female's survival chances to the next breeding season.
		In the current version of the model S_function is fixed for all females.
		
		Returns: 
			P_f : float
				value that represents the probability of survival
		"""

		investment = self.care + self.eggs
		P_f = np.exp(-1*(m_0 + m_r*investment**beta))
		return P_f
	def S_inseason(self):
		"""Calculates within-season survival probability of a female.
		If a female invests too much in care and egg production, it results in elevated mortality risk due to energy loss.
		S_inseason is used for fitness calculations and affects offspring survival.
		In the current version of the model S_inseason is fixed for all females.
		
		Returns: 
			S_f : float
				value that represents the probability of within-season survival 
		"""

		investment = self.care + self.eggs
		S_f = np.exp(-0.3*(m_0 + m_r*investment**beta))
		return S_f
	def set_offspring(self, partner_ID):
		self.offspring_list.append([self.fitness, self.epc, self.breeding_seasons, partner_ID])

class male():
	"""A class that represents male individuals 
	
	Attributes:
	nest_number: int
		position in a social network (social network == neighbourhood == tree)
		females and males occupying the same nest are referred to as a social pair
	tree_number: int
		identifies the social network number within the population
	global_nest_number: int
		global position in the population that consists of many social networks
	neighbours: list
		a list of neighbours that have connections with the individual
		indicates potential extra-pair mates
	f_epc: float
		extra-pair copulation rate (EPC) 
		evolving gene, not expressed in males, but present in their genome
	care_a: float
		slope of the care reaction norm (c_1)
		evolving gene, expressed only in males
	care_b: float
		intercept of the male care reaction norm (c_2); 
		evolving gene, expressed only in males
	defence_a: float
		slope in the X direction of the male defence reaction norm (d_1); 
		evolving gene, expressed only in males
	defence_b: float
		slope in the Y direction of the male defence reaction norm (d_2); 
		evolving gene, expressed only in males
	defence_c: float
		third coefficient of the male defence reaction norm (d_3); 
		evolving gene, expressed only in males
	border_defences: list of lists
		for each neighbour, a sub-list with neighbour's global_nest_number and partial territorial defence of the focal male is saved
	territory: float
		size of a territory secured by a male for his social partner and offspring in the nest
	mother: int
		global nest number of an individual's mother
	father: int
		global nest number of an individual's father
	offspring_list: list of lists 
		each sub-list contains paternity score, number of the breeding season and global_nest_number of the mother 
		the above values are appended separately for each offspring
	breeding_seasons : int
		age of the male, as a number of breeding seasons
	paternity_score : float
		probability of siring offspring, fitness equivalent for males
		calculated as a sum of mating probabilities with a social female (within-pair copulations) and all the neighbours (EPCs)
	"""

	def __init__(self): 
		self.care_a, self.care_b = 0., 0.5
		self.defence_a, self.defence_b, self.defence_c =  0., 0., 0.5
		self.f_epc = starting_epc
		self.mother = -1
		self.father = -1.
		self.offspring_list = []
		self.border_defences = []
		self.breeding_seasons = 1
		self.territory = 0.
		self.paternity_score = 0.

	def set_global_number(self, num): 
		self.global_nest_number = num
	def get_global_number(self):
		return self.global_nest_number
	def set_tree(self, num):
		self.tree_number = num
	def set_nest_number(self, num):
		self.nest_num = num
	def get_nest_number(self):
		return self.nest_num
	def set_neighbours(self, nlist):
		self.neighbours = []
		for n in nlist:
			self.neighbours.append(self.tree_number*individuals+n)
	def get_neighbours(self):
		return self.neighbours
	def parents(self, mother, father):
		"""Passes parameters to other functions that facilitate genetic inheritance.

			Parameters:
				mother : num
					global_nest_number of the mother
				father : num
					global_nest_number of the mother
		"""

		self.mother, self.father = mother, father
		self.set_care_params(female_parents[mother].get_male_care(), male_parents[father].get_care_params())
		self.set_defence_params(female_parents[mother].get_male_defence(), male_parents[father].get_defence_params())
		self.set_f_epc(female_parents[mother].get_epc(), male_parents[father].get_female_epc())
	def set_care_params(self, list1, list2):
		"""Calculates genetic coefficients for the reaction norm (c_1, c_2,).
		New values are drawn from a normal distribution with variance = 0.05.

		Parameters:
			list1 : list
				genetic values of c_1 and c_2 inherited from mother
			list2 : list
				genetic values of c_1 and c_2 inherited from father
		Returns: 
			care_a : float
				new slope of the male care reaction norm (c_1)
			care_b : float
				new intercept of the male care reaction norm (c_2)				
		"""

		par1 = population_epc_mean
		ang1 = math.atan(list1[0])
		ang2 = math.atan(list2[0])
		mean_a = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_a, 0.05)
		self.care_a = math.tan(new_angle)
		Y_parent1 = list1[0]*par1 + list1[1]
		Y_parent2 = list2[0]*par1 + list2[1]
		Y_offspring = (Y_parent1+Y_parent2)/2.
		new_b = Y_offspring - self.care_a*par1 
		self.care_b = np.random.normal(new_b, 0.05)
	def get_care_params(self):
		return [self.care_a, self.care_b]
	def set_defence_params(self, list1, list2):
		"""Calculates genetic coefficients for the reaction norm (d_1, d_2, d_3).
		New values are drawn from a normal distribution with variance = 0.05.

		Parameters:
			list1 : list
				genetic values of d_1, d_2 and d_3 inherited from mother
			list2 : list
				genetic values of d_1, d_2 and d_3 inherited from father
		Returns: 
			defence_a : float
				new slope in the X direction of the male defence reaction norm (d_1)
			defence_b : float
				new slope in the Y direction of the male defence reaction norm (d_2)
			defence_c : float
				new intercept of the male defence reaction norm (d_3)
		"""

		par1, par2 = population_epc_mean, (1.-population_epc_mean)/nmax
		ang1 = math.atan(list1[0])
		ang2 = math.atan(list2[0])
		mean_a = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_a, 0.05)
		self.defence_a = math.tan(new_angle)
		ang1 = math.atan(list1[1])
		ang2 = math.atan(list2[1])
		mean_b = (ang1+ang2)/2.
		new_angle = np.random.normal(mean_b, 0.05)
		self.defence_b = math.tan(new_angle)
		Z_val1 = list1[0]*par1 + list1[1]*par2 + list1[2]
		Z_val2 = list2[0]*par1 + list2[1]*par2 + list2[2]
		Z_val_offspring = (Z_val1+Z_val2)/2.
		new_c = Z_val_offspring - self.defence_a*par1 - self.defence_b*par2 
		self.defence_c = np.random.normal(new_c, 0.05)
	def get_defence_params(self):
		return [self.defence_a, self.defence_b, self.defence_c]
	def calculate_border_defence(self):
		"""Calculates partial territorial defence (delta) directed towards each neighbour.
			delta is calculated based on the social female's EPC value and EPC of the female neighbour.

		Returns: 
			a list of lists with partial territorial defences
			for each neighbour, a sub-list with their global_nest_number and partial territorial defence is saved
		"""

		for i in self.neighbours:
			arg1 = female_parents[i].get_epc()
			arg2 = 1 - female_parents[self.global_nest_number].get_epc()
			delta = self.defence_a*arg1 + self.defence_b*arg2 + self.defence_c
			if delta <= 0.:
				delta = 0.
			self.border_defences.append([i, delta])
	def get_expressed_d(self):
		"""Calculates D (expressed defence), the average of partial territorial defence

		Returns: 
			defence : float
		"""

		self.defence = 0.
		for d in range(0, len(self.border_defences)):
			self.defence += self.border_defences[d][1]
		self.defence = self.defence / len(self.neighbours)
		return self.defence
	def get_expressed_c(self):
		"""Calculates C (expressed care).Response to social female behaviour (EPC).

		Returns: 
			care : float
		"""

		self.care = self.care_a * female_parents[self.global_nest_number].get_epc() + self.care_b
		if self.care <= 0.:
			self.care = 0.
		return self.care
	def calculate_territory(self):
		"""Calculates territory secured by a social male based on his partial territorial defences (deltas) accross all the territory borders.

		Returns: 
			territory : float
		"""

		for i in self.border_defences:
			neighbour = i[0]
			delta_nm = i[1]
			for j in male_parents[neighbour].border_defences:
				if j[0] == self.global_nest_number:
					delta_mn = j[1]
					break
			t = (2./nmax) * delta_nm * self.S_inseason() / (delta_nm * self.S_inseason() + delta_mn * male_parents[neighbour].S_inseason() + floaters)
			self.territory += t
		for e in range(0, nmax-len(self.neighbours)):
			self.territory += 1./nmax
		return self.territory
	def S_function(self):
		"""Calculates between-seasons survival probability of a male.
		If a male invests too much in care and territorial defence, it results in elevated mortality risk due to energy loss.
		S_function is used to assess male's survival chances to the next breeding season.
		
		Returns: 
			P_m : float
				value that represents the probability of survival 
		"""

		investment = self.get_expressed_c() + self.get_expressed_d()
		P_m = np.exp(-1*(m_0 + m_r*investment**beta))
		return P_m
	def S_inseason(self):
		"""Calculates within-season survival probability of a male.
		If a male invests too much in care and territorial defence, it results in elevated mortality risk due to energy loss.
		S_inseason is used for female fitness calculations and affects offspring survival.
		
		Returns: 
			S_m : float
				value that represents the probability of within-season survival 
		"""

		investment = self.get_expressed_c() + self.get_expressed_d()
		self.S_m = np.exp(-0.3*(m_0 + m_r*investment**beta))
		return self.S_m
	def set_f_epc(self, val1, val2):
		"""Calculates genetic value for a female trait, x (EPC), from a normal distribution with variance = 0.05.
		If the new EPC value is below 0 or above 1, the function returns 0 or 1.

		Parameters:
			val1 : float
				genetic value of x inherited from mother
			val2 : float
				genetic value of x inherited from father
		Returns: 
			f_epc : float
		"""

		mean = (val1+val2)/2.
		self.f_epc = nrand.normal(mean, 0.05)
		if self.f_epc < 0.:
			self.f_epc = 0.
		if self.f_epc > 1.:
			self.f_epc = 1.
	def get_female_epc(self):
		return self.f_epc
	def calculate_paternity_score(self):
		self.paternity_score += (1. - female_parents[self.global_nest_number].get_epc())*female_parents[self.global_nest_number].get_fitness(minf, maxf)
		for n in self.neighbours:
			self.paternity_score += (female_parents[n].get_epc() / len(female_parents[n].neighbours))* female_parents[n].get_fitness(minf, maxf)
	def get_paternity_score(self):
		return self.paternity_score
	def set_offspring(self, partner_ID):
		self.offspring_list.append([self.get_paternity_score(), self.breeding_seasons, partner_ID])

class stats():
	""" Statistics
	Saving population parameters and stats. Trait evolution trajectories (epc, expressed care, expressed defence) 
		are saved after each breeding season (average and standard deviation).
	Detailed information about all genetic components (c_1, c_2, d_1, d_2, d_3), 
		for all individuals are saved every 100 generations.
	Optional: male and female performance.
	"""

	def __init__(self):
		self.mean_epc, self.sd_epc = [], []
		self.mean_fit, self.sd_fit = [], []
		self.mean_c, self.sd_c = [], []
		self.mean_d, self.sd_d = [], []
		self.carea, self.careb = [], []
		self.carea_f, self.careb_f = [], []
		self.defa, self.defb, self.defc = [], [], []
		self.defa_f, self.defb_f, self.defc_f = [], [], []
		self.epc_m, self.epc_f = [], []

	def input_males(self, table):
		self.males_list = table
		self.popul_D()
		self.popul_C()

	def input_females(self, table):
		self.females_list = table
		self.popul_EPC()
		self.popul_fitness()

	def input_details(self, tab1, tab2):
		self.females_list = tab1
		self.males_list = tab2
		self.care_params()
		self.defence_params()
		
	def popul_EPC(self):
		"""Calculates average EPC level expressed in the population.
		"""

		temp_table = []
		for i in self.females_list:
			temp_table.append(i.epc)	
		self.mean_epc.append(sta.mean(temp_table))
		self.sd_epc.append(sta.stdev(temp_table))

	def popul_fitness(self):
		"""Calculates average fitness value in the population.
		"""

		temp_table = []
		for i in self.females_list:
			temp_table.append(i.fitness)
		self.mean_fit.append(sta.mean(temp_table))
		self.sd_fit.append(sta.stdev(temp_table))

	def popul_D(self):
		"""Calculates average male defence expressed in the population.
		"""

		temp_table = []
		for i in self.males_list:
			temp_table.append(i.defence)
		self.mean_d.append(sta.mean(temp_table))
		self.sd_d.append(sta.stdev(temp_table))

	def popul_C(self):
		"""Calculates average male care expressed in the population.
		"""

		temp_table = []
		for i in self.males_list:
			temp_table.append(i.care)
		self.mean_c.append(sta.mean(temp_table))
		self.sd_c.append(sta.stdev(temp_table))

	def care_params(self):
		"""Extracts genetic values of c_1, c_2 and epc from males and females.
		"""

		for i in self.males_list:
			self.carea.append(i.get_care_params()[0])
			self.careb.append(i.get_care_params()[1])
			self.epc_m.append(i.get_female_epc())
		for i in self.females_list:
			self.carea_f.append(i.get_male_care()[0])
			self.careb_f.append(i.get_male_care()[1])
			self.epc_f.append(i.get_epc())			

	def defence_params(self):
		"""Extracts genetic values of d_1, d_2 and d_3 from males and females.
		"""
		for i in self.males_list:
			self.defa.append(i.get_defence_params()[0])
			self.defb.append(i.get_defence_params()[1])
			self.defc.append(i.get_defence_params()[2])
		for i in self.females_list:
			self.defa_f.append(i.get_male_defence()[0])
			self.defb_f.append(i.get_male_defence()[1])
			self.defc_f.append(i.get_male_defence()[2])

	def male_performance(self):
		for i in self.males_list:
			self.male_performance.append(i.offspring_list)

	def female_performance(self):
		for i in self.females_list:
			self.female_performance.append(i.offspring_list)

	def stats_out(self):
		"""After each breeding season, following features are saved:
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

		return [self.mean_epc[-1], self.sd_epc[-1], self.mean_c[-1], self.sd_c[-1], self.mean_d[-1], self.sd_d[-1]]

	def stats_details_out(self):
		"""Every 100 breeding seasons, following features are saved:
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
		
		num = trees*individuals
		return [self.carea[-num:], self.careb[-num:], self.carea_f[-num:], self.careb_f[-num:], \
			self.defa[-num:], self.defb[-num:], self.defc[-num:], self.defa_f[-num:], \
			self.defb_f[-num:], self.defc_f[-num:], self.epc_f[-num:], self.epc_m[-num:]]

# NETWORK TYPES AND OTHER FUNCTIONS
def small_world_network(nests, k, p):
	"""Watts-Strogatz (small world) network
	Uncomment lines 736-738 to plot the network

	Parameters:
		nests : number of nodes
		k : initial number of neighbours (k-1 if k is odd)
			k is the max number of neighbours after rewiring
		p : rewiring probability; set to zero in our case
	Returns:
		neighbours for all the nodes (as a list of lists) and k
	"""

	edges = []
	neighbours_list = [[] for x in range(0, nests)]
	connections = []
	G = nx.connected_watts_strogatz_graph(nests, k, p, tries=100, seed=None)
	for e in G.edges():
		neighbours_list[e[0]].append(e[1])
		neighbours_list[e[1]].append(e[0])
	#G = nx.watts_strogatz_graph(nests, k, p)
	# nx.draw(G, with_labels=True)
	# plt.show()
	count = []
	for n in neighbours_list:
		count.append(len(n))
	return neighbours_list, max(count)

def ladder_graph(nests):
	"""Ladder network 
	Creates a network with even number of nests, each with three connections
	Uncomment lines 760-762 to plot the network

	Parameters:
		nests : number of nodes/2
	Returns:
		neighbours for all the nodes (as a list of lists) and number of connections (3)
	"""

	neighbours_list = [[] for x in range(0, 2*nests)]
	G = nx.circular_ladder_graph(nests)
	for e in G.edges():
		neighbours_list[e[0]].append(e[1])
		neighbours_list[e[1]].append(e[0])
	# nx.draw(G, with_labels=True)
	# plt.plot(1, 1)
	# plt.show()
	return neighbours_list, 3

def complete(nests):
	"""Complete network
	Creates a fully connected network
	Uncomment lines 781-783 to plot the network

	Parameters:
		nests : number of nodes
	Returns:
		neighbours for all the nodes (as a list of lists) and number of connections
	"""

	neighbours_list = [[] for x in range(0, nests)]
	G = nx.complete_graph(nests)
	for e in G.edges():
		neighbours_list[e[0]].append(e[1])
		neighbours_list[e[1]].append(e[0])
	# nx.draw(G, with_labels=True)
	# plt.plot(1, 1)
	# plt.show()
	count = []
	for n in neighbours_list:
		count.append(len(n))
	return neighbours_list, max(count)

def floaters_inseason():
	investment = 1.8
	fl = 0.05 * np.exp(-0.3*(m_0 + m_r*investment**beta))
	return fl	

### Simulation parameters:
# mortality coefficient (m_0)								: 0.4	OR	0.2	
# female reproductive investment (egg production + care)	: 1.2 (0.6 + 0.6)	OR	1.8 (0.9 + 0.9)
# length of the simulation (years)							: 160 000	OR	100 000
# starting EPC value										: 0.0	OR	0.2	OR	0.4	OR	0.6

# Declare number of years per run, number of networks (trees) and individuals of each sex per tree
years = 100000
trees = 220
individuals = 24
m_0 = 2.0
m_r = 0.1 
starting_epc = 0.0
female_care, eggs = 0.9, 0.9
fname = 'p1_s1_'
beta = 3.
alpha = 0.7
gamma = 0.7

if __name__ == "__main__":
	start_time = time.time()
	female_parents, male_parents = [], []
	non_breeders_f, non_breeders_m = [], []
	female_offspring, male_offspring = [], []
	lucky_females, lucky_males = [], []
	floaters = floaters_inseason()
	
	# 1. START: repeat simulation for a given number of years
	for y in range(0, years):
		if y%100 == 0:
			print(y)
		if y%20000== 0:
			output_data = stats()
		find_fitness = []
		population_epc_mean = 0.
		fitness_sum = 0

		# 2a. Year zero: create lists of female parents and male parents
		if y == 0:
			for i in range(0, trees*individuals):
				female_parents.append(female())
				male_parents.append(male())

		# 2b. Create new lists consisting of males and females that enter the new breeding season
		# The lists contain reshuffled offspring, and "lucky males" and "lucky females" (survivors from the previous seasons)
		# "lucky males" and "lucky females" have priority when establishing a nest
		# The offspring that can't establish a nest (when all the nests are already occupied) are called "non-breeders"
		else:
			female_parents, male_parents = [], []
			new_male_parents, new_female_parents = trees*individuals - len(lucky_males), trees*individuals - len(lucky_females)
			female_parents = lucky_females
			male_parents = lucky_males
			nrand.shuffle(female_offspring)
			nrand.shuffle(male_offspring)
			male_parents += male_offspring[0:new_male_parents]
			female_parents += female_offspring[0:new_female_parents]
			non_breeders_f, non_breeders_m = female_offspring[new_female_parents:], male_offspring[new_male_parents:]
			lucky_males, lucky_females = [], []
			female_offspring, male_offspring = [], []
			nrand.shuffle(female_parents) 
			nrand.shuffle(male_parents) 

		# 3. Assign global nest numbers to each individual
		for i in range(0, len(female_parents)):
			female_parents[i].set_global_number(i)
			male_parents[i].set_global_number(i)

		# 4. For each individual: find neighbours and tree number 
		for t in range(0, trees):
			neighbours, nmax = ladder_graph(int(individuals/2))
			# alternative networks:
			# neighbours, nmax = small_world_network(individuals, 6, 0.)
			# neighbours, nmax = complete(individuals)
			for n in range(0, len(neighbours)):
				female_parents[individuals*t+n].set_tree(t)
				male_parents[individuals*t+n].set_tree(t)
				female_parents[individuals*t+n].set_neighbours(neighbours[n])
				male_parents[individuals*t+n].set_neighbours(neighbours[n])

		# 5. For each male: calculate partial territorial defence aimed towards all neighbours, and within-season survival
		for m in male_parents:
			m.calculate_border_defence()
			m.S_inseason()

		# 6. Based on male reproductive investment (C and D), calculate and normalize female fitness
		# Calculate average EPC level expressed in the whole population
		for f in female_parents:
			find_fitness.append(f.W_function())
		minf, maxf = min(find_fitness), max(find_fitness)
		for f in female_parents:
			fitness_sum += f.get_fitness(minf, maxf)
			population_epc_mean += f.get_epc()
		population_epc_mean = population_epc_mean/len(female_parents)

		# 7. Among breeders and non-breeders:
		# find individuals that survive to the next breeding season - "lucky males" and "lucky females"
		# This has to be done prior to breeding, number of surviving individuals defines the number of offspring 
		for m in male_parents:
			m.calculate_paternity_score()
			random_number = nrand.uniform(0., 1.)
			if random_number <= m.S_function():
				lucky_males.append(m)

		for f in female_parents:
			random_number = nrand.uniform(0., 1.)
			if random_number <= f.S_function():
				lucky_females.append(f)

		for n in non_breeders_m:
			random_number = nrand.uniform(0., 1.)
			S_function = np.exp(-1*(m_0 + 1)) # not breeding, reproductive investment == 0
			if random_number <= S_function:
				lucky_males.append(n)

		for k in non_breeders_f:
			random_number = nrand.uniform(0., 1.)
			S_function = np.exp(-1*(m_0 + 1)) # not breeding, reproductive investment == 0
			if random_number <= S_function:
				lucky_females.append(k)

		# 8. Find the rare sex and calculate offspring number (optimization)
		if len(lucky_females) > len(lucky_males):
			rare_sex = len(lucky_males)
		else:
			rare_sex = len(lucky_females)

		offspring_total = trees*individuals - rare_sex

		# 9. Calculate cumulative fitness of all females in the tree
		for l in range(0, trees):
			control_sum = 0
			for f in female_parents[l*individuals:l*individuals+individuals]:
				if control_sum == 0:
					cumulative_fitness = np.array(f.get_fitness(minf, maxf))
					control_sum+=f.get_fitness(minf, maxf)
				else:
					cumulative_fitness = np.append(cumulative_fitness, f.get_fitness(minf, maxf))
					control_sum+=f.get_fitness(minf, maxf)
			cumulative_fitness = np.cumsum(cumulative_fitness)/control_sum

			# 10. Calculate the number of male- and female-offspring per tree, based on the sum of female fitness
			n_offspring = int(round((offspring_total*control_sum/fitness_sum),0)) + 1 # rounding rule

			# 11. Among all females in a tree: choose the female who mates based on her fitness measure
			# female_ID = global_nest_number of a chosen female
			for k in range(0, n_offspring*2):
				if cumulative_fitness[-1] > 1.:
					rand = 1.
				else:
					rand = nrand.uniform(0, cumulative_fitness[-1])
				for i in range(0, len(cumulative_fitness)):
					if rand in portion.closedopen(0, cumulative_fitness[i]):
						female_ID = l*individuals+i 
						break
				# 12. Find a male candidate based on female's epc value and her number of neighbours
				# male_ID = global_nest_number of a sire
				if k < n_offspring:
					# 13a. Create object 'male' and add it to the male_offspring list
					male_offspring.append(male())
					male_ID = female_parents[female_ID].choose_mate()
					male_offspring[-1].parents(female_ID, male_ID)
					female_parents[female_ID].set_offspring(male_ID)
					male_parents[male_ID].set_offspring(female_ID)
				else:
					# 13b. Create object 'female' and add it to the female_offspring list
					female_offspring.append(female())
					male_ID = female_parents[female_ID].choose_mate()
					female_offspring[-1].parents(female_ID, male_ID)
					female_parents[female_ID].set_offspring(male_ID)
					male_parents[male_ID].set_offspring(female_ID)

		# 14. Statistics
		if y < 20000:
			file_stats = open(str(fname)+"0_traits.pickle", "ab")
			details = open(str(fname)+"0_details.pickle", "ab")
		elif y >= 20000 and y < 40000:
			file_stats = open(str(fname)+"1_traits.pickle", "ab")
			details = open(str(fname)+"1_details.pickle", "ab")
		elif y >= 40000 and y < 60000:
			file_stats = open(str(fname)+"2_traits.pickle", "ab")
			details = open(str(fname)+"2_details.pickle", "ab")
		elif y >= 60000 and y < 80000:
			file_stats = open(str(fname)+"3_traits.pickle", "ab")
			details = open(str(fname)+"3_details.pickle", "ab")
		elif y >= 80000 and y < 100000:
			file_stats = open(str(fname)+"4_traits.pickle", "ab")
			details = open(str(fname)+"4_details.pickle", "ab")
		elif y >= 100000 and y < 120000:
			file_stats = open(str(fname)+"5_traits.pickle", "ab")
			details = open(str(fname)+"5_details.pickle", "ab")
		elif y >= 120000 and y < 140000:
			file_stats = open(str(fname)+"6_traits.pickle", "ab")
			details = open(str(fname)+"6_details.pickle", "ab")
		elif y >= 140000 and y < 160000:
			file_stats = open(str(fname)+"7_traits.pickle", "ab")
			details = open(str(fname)+"7_details.pickle", "ab")
		elif y >= 160000 and y < 180000:
			file_stats = open(str(fname)+"8_traits.pickle", "ab")
			details = open(str(fname)+"8_details.pickle", "ab")
		else:
			file_stats = open(str(fname)+"9_traits.pickle", "ab")
			details = open(str(fname)+"9_details.pickle", "ab")

		# Each year: save general stats
		output_data.input_males(male_parents)
		output_data.input_females(female_parents)
		pickle.dump(output_data.stats_out(), file_stats)
		file_stats.close()
		# Every 100 years: save detailed stats
		if y%100 == 0:
			output_data.input_details(female_parents, male_parents)
			pickle.dump(output_data.stats_details_out(), details)
			details.close()

		# 15. Update/ restore default values for individuals surviving to the next breeding season
		for lm in lucky_males:
			lm.breeding_seasons += 1
			lm.border_defences = []
			lm.territory = 0.
			lm.paternity_score = 0.

		for lf in lucky_females:
			lf.breeding_seasons += 1
			lf.territory = 0.
			lf.fitness = 0.

	# 16. End; Print simulation time
	print('time: '+str((time.time() - start_time)))