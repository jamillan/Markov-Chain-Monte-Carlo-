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
	def __init__(self , filepath = None, seed=1, totalsteps=0,initnode=0,finalnode=0 , timeskips =10000):
		#Initalize Instance of Markovian Chain
		self.MATRIX_SET =False;
		if filepath != None:
		  self.adjMatrix = np.genfromtxt(filepath)
		  self.MATRIX_SET= False ; 

		self.rand =random.Random()
		self.rand.seed(seed)
		self.timeskips = int(timeskips)
		self.totalsteps = int(totalsteps)
		self.initnode = int(initnode)
		self.finalnode = int(finalnode)
		self.globproperty = 0.0

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
				print("\t| steps remaining : " , steps , "/" , self.totalsteps  ,"\t| time remaining : " , (self.totalsteps - steps) * ETS , "seconds") 
				loop = loop + 1
				freq = 0				

			#probability for trial move
			cutoff = self.rand.random();
			prob = 0;
			#print(self.adjMatrix) 
			#quit() 
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

		print("****Simulation Complete****")

		return counts/self.totalsteps

	def diffussion(self , initstate, finalstate):
	# run Markovian Chain for a targeted node distribution
	# variables for probability calculations
		counts = 0.
		freq = 0.
		loop = 1
		diffsteps = 0.
		listDiff = []	
    

		FinalNode = int(initstate)


		if type(initstate) is not int: 
			print("Warning : the initial node argument was not integer but it we try to convertwas convert integer") 
			initstate = int(initstate)


	 #Check if  initial and final nodes are contained in the transition probability  matrix
		if((finalstate+1) > len(self.adjMatrix)):
			raise RuntimeError("ERROR: Final node label in larger that size of adjacent Matrix: " + str(len(self.adjMatrix)))

		if((initstate + 1)  > len(self.adjMatrix)):
			raise RuntimeError("ERROR: Final node label in larger that size of adjacent Matrix: " + str(len(self.adjMatrix)))
	 #Check if starting from final Node and throw error if so !
		if(initstate == finalstate):
			raise RuntimeError("ERROR: Starting from final targeted final node")
		# Loop for trial moves 

		for steps in range(self.totalsteps):


			if freq == 0. :
				t_init = time.time()

			loop = 1
			if steps % self.timeskips == 0 and steps > 0:
				ETS = (time.time() - t_init) / self.timeskips
				print("\t| steps remaining : " , steps , "/" , self.totalsteps  ,"\t| time remaining : " , (self.totalsteps - steps) * ETS , "seconds" ) 
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
					diffsteps += 1.0;

					if type(finalstate)  is int :
					  if FinalNode  == finalstate :
						  listDiff.append(diffsteps)

					if type(finalstate)  is list:
					  if FinalNode  in finalstate :
						  listDiff.append(diffsteps)
					

					#print "State | " + str(FinalNode)
					break

			if type(finalstate) is int: 
			  if FinalNode == finalstate: 
				  FinalNode = initstate;
				  diffsteps = 0.0;

			if type(finalstate) is list: 
			  if FinalNode in finalstate: 
				  FinalNode = initstate;
				  diffsteps = 0.0;
		

		#return probability of targeted state

		return listDiff

	def setinitnode(self,init):
		#set initial node
		if init < (len(self.adjMatrix[0])-1): 
			self.initnode = init

		else:
			print("wrong nbr of states: total nbr " + str(len(self.adjMatrix[0])-1))

	def setfinalnode(self, final):
	#ser final node
		if init < (len(self.adjMatrix[0])-1): 
			self.initnode = final

		else:
			print("wrong nbr of states: total nbr " + str(len(self.adjMatrix[0])-1))
	
	def raisematrix(self, nth):
		#raise transition matrix to the "nth" power
		powmatrix = LA.matrix_power(self.adjMatrix, nth)
		return powmatrix;
			
	def set_adjmatrix(self,A):
		self.adjMatrix = np.array(A,dtype=float)
		self.MATRIX_SET = True ;
		print(self.adjMatrix) 

# compute transition function
	def compute_transition(self,node):
		if self.MATRIX_SET ==False:
		  raise RuntimeError("Forgot to set transition Matrix")

		if self.MATRIX_SET ==True:
		  cutoff = self.rand.random();
		  prob = 0;

		if type(node) is not int:
		  node = int(node);
		  print("Node argument is not of integer type: Attemp to transform to integer was performed but it may have failed")

		for j in range(len(self.adjMatrix)):
		  prob+= self.adjMatrix[node,j];
		  if prob >= cutoff:
		   return j;


		
	def sampling_property(self,A):
		self.property = A
		if len(A) != len(self.adjMatrix):raise RuntimeError("sampling should be of size " + str(len(self.adjMatrix))); 


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
	runner = MC(filepath = var[0], seed = var[1],totalsteps = var[2], initnode = var[3], finalnode = var[4], timeskips = var[5])
	#get equil_dist
	equil_dist  = runner.run()
	return equil_dist , runner

if __name__ == '__main__':
	dir = os.getcwd()
	

	# Example of Random Walker 
	# Read Parameters
	totalsteps = 1e6    #Total Steps
	seed =  3							 # Seed(s) for Random number(s)
	skipsteps = 10000			 # Profiling time
	node_0  =  0  				 # Initial Node
	node_1  = 1   				 # List of nodes to obtain its equilibrium Distributions
	


  # Instance of Markovian Chain
	# File randwalk.dat contains transition matrix
 	#init node 0
	#let's get equil distribution of node 1
	init_node = 0
	equil_node = 1
	r_seed = 4
	randwalker = MC(  seed = r_seed , totalsteps = totalsteps ,initnode = init_node , finalnode = equil_node,timeskips =  10000)
	matrix = np.array([[0.3,0.7],[0.5,0.5]])
	randwalker.set_adjmatrix(matrix)
	j = randwalker.compute_transition(0)
	print(j)
	quit()

	#raise transition matrix to the 100th power
	A = randwalker.raisematrix(100)

	print("equil distribution matrix" , A)

	# Let's do some sampling
	# Lets say all nodes have age = 1
	prop = np.ones(2)
	#but node 0 is 10 years old in reality
	prop[0]=10.0
	
	#pass sampling property to the MC
	randwalker.sampling_property(prop)
  
	#Sanple
	sampling =True
	prob_node = randwalker.run(sampling  = sampling)

	#Print the Global Age
	print("Global Age : ",randwalker.globalproperty())
	print( "equil distribution of node 1: " , prob_node)

	#Get average number steps from node 1 to node 2
	diff_1_2 = randwalker.diffussion(int(0),int(1))
	print( "average steps from node 1 to node 2 : " , np.mean(diff_1_2))


  # Lets count average number of edges aunt visits in a reversible walk between opposite vertex vertex 
  # Imagine an aunt starts at some vertex (node = 0) on a cube and randomly travels to opposite vertex (node =6)  
  # How many vertices visits in average before completing the full trip?

	matrix = np.genfromtxt("aunt_edges.dat",dtype=float)
	randwalker.set_adjmatrix(matrix) 
	j =randwalker.compute_transition(0) 
	diff_1_2 = randwalker.diffussion(int(0),int(6))
	print( "Aunt visits an average number of edges equivalent to: " , np.mean(diff_1_2))



  #parallel Simulations
	#Create list of target nodes and other parameters

	var=[[ dir + "/matrix.dat" , 0 , totalsteps ,0 , 0 , skipsteps]]
	var1=[ dir + "/matrix.dat" , 1 , totalsteps ,1 ,0 , skipsteps]

	var.append(var1)

	#example of  how to run in parallel
	#USe encapsulate function to run in parallel
	matrix = np.array([[0.3,0.7],[0.5,0.5]])
	pool = Pool(processes = multiprocessing.cpu_count())
	runners =  pool.map(encapsulate,var)


	for i in range(len(runners)):
		
		print("printing distribution node " + str(i) + " = " ,runners[i][0])
	
	





