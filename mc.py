from __future__ import division 
import sys
import numpy as np
from numpy import linalg as LA
import pandas as pd
import random
import multiprocessing
from multiprocessing import Pool
import time 
import json
import os



class MC:
	def __init__(self , filepath, seed, totalsteps,initnode,finalnode , timeskips =10000):
		#Initalize Instance of Markovian Chain
		self.adjMatrix = np.genfromtxt(filepath)
		self.rand =random.Random()
		self.rand.seed(seed)
		self.timeskips = timeskips
		self.totalsteps = totalsteps
		self.initnode = initnode
		self.finalnode = finalnode
		self.globproperty = 0

		#sanity check 

		

	def read(filepath):
		self.adjMatrix = np.genfromtxt(filepath)

	def run(self,sampling=False):
	# run Markovian Chain for a targeted node distribution
	# variables for frequency calculations
		counts = 0.
		freq = 0.
		FinalNode = self.initnode
		loop = 1
		globalproperty =0
		for steps in range(self.totalsteps):
			#print "steps | " + str(steps)
			if freq == 0. :
				t_init = time.time()

			loop = 1
			if steps % self.timeskips == 0 and self.totalsteps > 0:
				ETS = (time.time() - t_init) / self.timeskips
				print "\t| steps remaining : " , steps , "/" , self.totalsteps  ,"\t| time remaining : " , (self.totalsteps - steps) * ETS , "seconds"  
				loop = loop + 1
				freq = 0				

			#probability for trial move
			cutoff = self.rand.random();
			prob = 0;
			#Trial Move
			for j in range(len(self.adjMatrix)):
				prob += self.adjMatrix[FinalNode,j];
				
				if prob >= cutoff:
					FinalNode = j

					#if sampling
					if sampling:
						globalproperty = globalproperty + self.property[j]

					if FinalNode == self.finalnode :
					# Calculate probability for targeted node
						counts= counts + 1.0
						freq = freq + 1.0

					#print "State | " + str(FinalNode)
					break
		
		#self.finalnode = FinalNode;

		#return probability of targeted state

		if sampling:
			self.globproperty = globalproperty / self.totalsteps

		print "****Simulation Complete****"

		return counts/self.totalsteps

	def diffussion(self , initstate, finalstate):
	# run Markovian Chain for a targeted node distribution
	# variables for probability calculations
		counts = 0.
		freq = 0.
		FinalNode = self.initnode
		loop = 1
		edges = 0. 
		diffsteps = 0.
		listDiff = []	
		# Loop for trial moves 

		for steps in range(self.totalsteps):

			#print "steps | " + str(steps)
			if(FinalNode == initstate): 
				edges  = 0
				diffsteps = 1

			if freq == 0. :
				t_init = time.time()

			loop = 1
			if steps % self.timeskips == 0 and steps > 0:
				ETS = (time.time() - t_init) / self.timeskips
				print  "\t| steps remaining : " , steps , "/" , self.totalsteps  ,"\t| time remaining : " , (self.totalsteps - steps) * ETS , "seconds"  
				loop = loop + 1
				freq = 0				

			#probability for trial move
			cutoff = self.rand.random();
			prob = 0;
			#Trail Move
			for j in range(len(self.adjMatrix)):
				prob += self.adjMatrix[FinalNode,j];
				
				if prob >= cutoff:
					FinalNode = j

					if FinalNode  == finalstate :
						diff = edges ; 
						listDiff.append(diff)
					
					else:
						edges = edges + 1.0
						diffsteps = diffsteps + 1.0


					if FinalNode == self.finalnode :
					# Calculate probability for targeted node
						counts= counts + 1.0
						freq = freq + 1.0

					#print "State | " + str(FinalNode)
					break
		

		#return probability of targeted state

		return listDiff

	def setinitnode(self,init):
		#set initial node
		if init < (len(self.adjMatrix[0])-1): 
			self.initnode = init

		else:
			print "wrong nbr of states: total nbr " + str(len(self.adjMatrix[0])-1)

	def setfinalnode(self, final):
	#ser final node
		if init < (len(self.adjMatrix[0])-1): 
			self.initnode = final

		else:
			print "wrong nbr of states: total nbr " + str(len(self.adjMatrix[0])-1)
	
	def raisematrix(self, nth):
		#raise transition matrix to the "nth" power
		powmatrix = LA.matrix_power(self.adjMatrix, nth)
		return powmatrix;
			
	def adjmatrix(self,A):
		self.adjMatrix = A
		
	def sampling_property(self,A):
		self.property = A 

	def globalproperty(self):
		return self.globproperty 

	def eigen(self):
		powmatrix = LA.matrix_power(self.adjMatrix, 1)
		powmatrix = np.transpose(powmatrix)
		eigenval,eigenvec = LA.eig(powmatrix)
		return eigenval,eigenvec

#wrapper function for parallel execution

def encapsulate(var,sampling=False):
	#create its own MC runner	
	runner = MC(var[0],var[1], var[2], var[3], var[4])
	#get equil_dist
	equil_dist  = runner.run()
	return equil_dist , runner

if __name__ == '__main__':
	dir = os.getcwd()
	

	# Example of Random Walker 
	# Read Parameters
	totalsteps = 10000    #Total Steps
	seed =  2							 # Seed(s) for Random number(s)
	skipsteps = 1000			 # Profiling time
	node_0  =  0  				 # Initial Node
	node_1  = 1   				 # List of nodes to obtain its equilibrium Distributions
	


  # Instance of Markovian Chain
	# File randwalk.dat contains transition matrix
 	#init node 0
	#let's get equil distribution of node 1
	init_node = 0
	equil_node = 1
	r_seed = 212
	randwalker = MC( r_seed , totalsteps ,init_node , equil_node, 10000)
	matrix = np.array([[0.5,0.5],[0.5,0.5]])
	randwalker.adjmatrix(matrix)

	#raise transition matrix to the 100th power
	A = randwalker.raisematrix(100)

	print "equil distribution matrix" , A

	# Let's do some sampling
	# Lets say all nodes have age = 1
	prop = np.ones(5)
	#but node 0 is 10 years old in reality
	prop[0]=10.0
	
	#pass sampling property to the MC
	randwalker.sampling_property(prop)
  
	#Sanple
	sampling =True
	prob_node = randwalker.run(sampling)

	#Print the Global Age
	print "Global Age : ",randwalker.globalproperty()
	print "equil distribution of node 1: " , prob_node

	#Get average number steps from node 1 to node 2
	diff_1_2 = randwalker.diffussion(1,2)
	print "average steps from node 1 to node 2 : " , np.mean(diff_1_2)

  #parallel Simulations
	#Create list of target nodes and other parameters
	var=[[ 0 , totalsteps ,1 , 0 , skipsteps]]
	var1=[ 1 , totalsteps ,2 ,1 , skipsteps]
	var.append(var1)

	#example of  how to run in parallel
	#USe encapsulate function to run in parallel
	pool = Pool(processes = multiprocessing.cpu_count())
	runners =  pool.map(encapsulate,var)

	for i in range(len(runners)):
		
		print "printing distribution node " + str(i) + " = " ,runners[i][0]
	
	





