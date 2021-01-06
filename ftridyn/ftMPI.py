#!/usr/bin/env python

from mpi4py import MPI
#import generate_ftridyn_input_gitr
#import analyze_ftridyn_simulations_gitr
#import ftridyn
import subprocess
import os
import time
#import subprocess
import numpy as np
import netCDF4
import getopt, sys
import pickle
from io import StringIO
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
WORKTAG=0
DIETAG=1

class species:
	def __init__(self, ZZ, M, BE, ED, EF, QU, DNS0, CK, E0, ALPHA0, QUBEAM, QUMAX, SBV):
		self.ZZ = ZZ					#Atomic number/charge
		self.M = M						#Mass in [u]
		self.BE = BE					#Binding energy in [eV]
		self.ED = ED					#Displacement energy in [eV]
		self.EF = EF					#Cutoff energy in [eV]
		self.QU = QU					#Concentration of species in homogeneous target
		self.DNS0 = DNS0			#Number density in [at/A^3]
		self.CK = CK					#Electronic stopping correction factor (0.0<CK<1.0)
		self.E0 = E0					#Incident energy in [eV]
		self.ALPHA0 = ALPHA0	#Incident angle in [deg] (0.0<ALPHA0<90.0)
		self.QUBEAM = QUBEAM	#Concentration in beam (0.0<QUBEAM<1.0)
		self.QUMAX = QUMAX		#Maximum allowed concentration in target (0.0<QUMAX<1.0)
		self.SBV = SBV				#Surface binding energy in [eV]
	#end def __init__
#end class species

class species_lookup_table:
	def __init__(self,dataFile):
		#cwd = os.path.dirname(os.path.realpath(__file__))
		print('dataFile',dataFile)
		self.table1 = open(dataFile,'r').readlines()[11:130]
		self.species_name = []
		self.ZZ = []
		self.M = []
		self.BE = []
		self.ED = []
		self.EF = []
		self.DNS0 = []
		self.CK = []
		self.SBV = []
		for i in range(len(self.table1)):
			line = self.table1[i]
			self.species_name.append(line[0:5].strip())
			self.ZZ.append(float(line[6:10].strip()))
			self.M.append(float(line[11:22].strip()))
			self.BE.append(float(line[114:120].strip()))
			self.ED.append(float(line[53:59].strip()))
			self.EF.append(float(line[62:68].strip()))
			self.DNS0.append(float(line[32:41].strip()))
			self.SBV.append(float(line[43:49].strip()))
		#end for
	#end def __init__

	def find_species(self,name,QU=0.0,CK=1.0,E0=0.0,ALPHA0=0.0,QUBEAM=0.0,QUMAX=0.0):
		species_index = self.species_name.index(name)
		return species(ZZ=self.ZZ[species_index],
			M=self.M[species_index], BE=self.BE[species_index],
			ED=self.ED[species_index], EF=self.EF[species_index],
			QU=QU, DNS0=self.DNS0[species_index], CK=CK, E0=E0, ALPHA0=ALPHA0,
			QUBEAM=QUBEAM, QUMAX=QUMAX, SBV=self.SBV[species_index])
	#end def find_species
#end def species_lookup

class sim_params:
	#simulation will hold a list of species and simulation parameters to be
	#declared in main().
	def __init__(self,name,sim_num,IFOUT,NH,IDOUT,IQOUT,NCP,IDREL,IQ0,IRC0,IRAND,
		JSP1,JSP2,JFRP,JNRM,FLC,INEL,IWC,IDIFF,CSBE,ANGMN,ANGMX,TT,TTDYN,NQX,DSF,
		IQXN,IQXX,IMCP,surf_name,species_list,output_pka,output_ska,output_prj):
		self.name = name 				#4-character file name, e.g., He_W, C_C_, NaNa, etc.
		self.sim_num = number_to_4_string(sim_num) 	#4-number simulation number, e.g., 0001, 0002
		self.IFOUT = int(IFOUT) 			#Number of projectiles after which info is written to terminal
		self.NH = int(NH)							#Number of incident histories (i.e. projectiles)
		self.IDOUT = int(IDOUT) 			#Integral data ouput after every IDOUT'th history
		self.IQOUT = int(IQOUT) 			#Additional .PR## profile file output after every IQOUT'th history
		self.NCP = int(NCP) 					#Number of components (i.e. species)
		self.IDREL = int(IDREL) 			#IDREL>0:Suppression of relaxation, IDREL<0:Suppression of relaxation and cascades
		self.IQ0 = int(IQ0) 					#IQ0=0:homogeneous initial composition, IQ0!=0:inhomogeneous initial composition from .LAY
		self.IRC0 = int(IRC0) 				#IRC0<0: Subthreshold atoms free, IRC0>=0: Subthreshold atoms bound
		self.IRAND = int(IRAND) 			#Random seed
		self.JSP1 = int(JSP1) 				#Suppression of recoils of type JSP1,...,JSP2
		self.JSP2 = int(JSP2) 				#JSP1=0: JSP2 ignored, all recoils traced
		self.JFRP = int(JFRP) 				#Frenkel Pair generation for JFRP,...,NCP
		self.JNRM = int(JNRM) 				#Normalzation of print output
		self.FLC = FLC 								#Fluence [1e16/cm2]
		self.INEL = int(INEL) 				#INEL=[1,2,3]:[inelastic nonlocal, local, equipartition]
		self.IWC = int(IWC) 					#Max order of weak projectile-target collisions {1,2,3}
		self.IDIFF = int(IDIFF) 			#IDIFF=[0,1]:[reemission,"diffusion"]
		self.CSBE = CSBE 				#Surface binding energy model (see documentation)
		self.ANGMN = ANGMN 			#Minimum angle for sputtered energy/angle distributions
		self.ANGMX = ANGMX 			#Maximum angle for sputtered energy/angle distributions
		self.TT = TT 						#Depth of surface
		self.TTDYN = TTDYN 			#Depth of dynamic surface (>=TT)
		self.NQX = NQX 					#Number of composition layers within F-TRIDYN
		self.DSF = DSF 					#Averaging depth for surface composition
		self.IQXN = IQXN 				#Limiting depth intervals for .PR## profile output
		self.IQXX = IQXX				#
		self.IMCP = int(IMCP) 				#Component for which moments should be calcluated

		self.surf_name = surf_name 				#Name of surface file
		self.species_list = species_list 	#List of species objects

		self.output_pka = output_pka #PKA output from F-TRIDYN
		self.output_ska = output_ska #SKA output from F-TRIDYN
		self.output_prj = output_prj #Projectile output from F-TRIDYN
		#todo: error checking, NCP=length(species_list), all parameters within bounds, etc.
	#end def __init__

	def print_input_file(self):
		input_file = open(self.name+self.sim_num+'.IN','w+')
		print(self.name+self.sim_num+'.IN', file=input_file)
		print(self.IFOUT,file=input_file)
		print(self.NH, self.IDOUT, self.IQOUT, self.NCP, self.IDREL, self.IQ0,
			self.IRC0, self.IRAND, self.JSP1, self.JSP2, self.JFRP, self.JNRM,
			file=input_file)
		print(self.FLC, self.INEL, self.IWC, self.IDIFF, self.CSBE, self.ANGMN,
			self.ANGMX, file=input_file)
		print(self.TT, self.TTDYN, self.NQX, self.DSF, self.IQXN, self.IQXX,
			self.IMCP, file=input_file)
		print(self.surf_name,self.output_pka, self.output_ska, self.output_prj, file=input_file)
		for i in range(self.NCP):
				current_species = self.species_list[i]

				print(current_species.ZZ, current_species.M, current_species.BE,
					current_species.ED, current_species.EF, current_species.QU,
					current_species.DNS0, current_species.CK, file=input_file)

				print(current_species.E0, current_species.ALPHA0,
					current_species.QUBEAM, current_species.QUMAX, end='',
					file=input_file)

				for j in range(self.NCP):
					print(' '+str(current_species.SBV), end='', file=input_file)
				#end for

				print(file=input_file)
		#end for
	#end def print_input_file
#end class simulation

class fractal_surface:
	def __init__(self, FD=1.0, width=400.0, iterations=1):
		self.FD = FD
		self.width = width
		if FD==1.0:
			self.iterations = 1
		else:
			self.iterations = iterations
		self.fgen()
		self.name = str(1)+'p'+str(np.int((FD-0.99999)*1000.0))
	#end def __init__

	def fgen(self,shape=[0,1,0,-1,0,-1,0,1,0]):
		b = np.sum(np.abs(shape))
		a = np.size(shape)-b
		beta = np.arccos((-a+(a+b)**(1.0/self.FD))/b)
		L = 1.0/(a+b*np.cos(beta))
		num_gen_seg = np.size(shape)
		self.num_gen_points = num_gen_seg+1

		self.x = np.zeros(self.num_gen_points)
		self.y = np.zeros(self.num_gen_points)
		segX = np.zeros(self.num_gen_points)
		segY = np.zeros(self.num_gen_points)
		self.gx = np.zeros(self.num_gen_points)
		self.gy = np.zeros(self.num_gen_points)

		self.x[0] = 0.0
		self.y[0] = 0.0

		for i in range(1,self.num_gen_points):
			self.x[i] = self.x[i-1] + L * np.cos(beta*shape[i-1])
			self.y[i] = self.y[i-1] + L * np.sin(beta*shape[i-1])
		#end for

		self.gx = self.x
		self.gy = self.y

		for i in range(2,self.iterations+1):
			counter = -1
			storeX = self.x
			storeY = self.y
			self.x = np.zeros(num_gen_seg**(i)+1)
			self.y = np.zeros(num_gen_seg**(i)+1)

			for j in range(np.size(storeX)-1):
				segX[0] = storeX[j]
				segY[0] = storeY[j]
				gamma = np.arctan2((storeY[j+1]-segY[0]),(storeX[j+1]-segX[0]))

				for k in range(1,num_gen_seg+1):
					coordX = L**i*np.cos(beta*shape[k-1])
					coordY = L**i*np.sin(beta*shape[k-1])
					segX[k] = segX[k-1] + coordX*np.cos(gamma) - coordY*np.sin(gamma)
					segY[k] = segY[k-1] + coordX*np.sin(gamma) + coordY*np.cos(gamma)
				#end for k

				if(j<np.size(storeX)-2):
					segTrim = -1
				else:
					segTrim = 0
				#end if

				for n in range(0,np.size(segX)+segTrim):
					counter = counter + 1
					self.x[counter] = segX[n]
					self.y[counter] = segY[n]
				#end for n
			#end for j
		#end for i
		self.num_surf_points = np.size(self.x)
		self.x,self.y = self.width*self.x,self.width*self.y #scale the final surface
	#end def fgen

	def print_fractal_surface_file(self):
		surface_file = open(self.name+'.surf','w+')
		#Surface file format is as follows:
		#number of points in generator, number of iterations, number of points in surface
		#fractal dimension
		#generator x points, generator y points
		#...
		#fractal x points, fractal y points
		#...
		print(self.num_gen_points, self.iterations, self.num_surf_points, file=surface_file)
		print(self.FD,file=surface_file)

		for i in range(self.num_gen_points):
			print(self.gx[i], self.gy[i], file=surface_file)
		#end for

		for i in range(self.num_surf_points):
			print(self.x[i], self.y[i], file=surface_file)
		#end for

	#end def print_fractal_surface_file
#end class fractal_surface

def number_to_4_string(number):
	string_number = str(number)
	#todo: try except blocks to capture error
	length_string = len(string_number)
	for i in range(4-length_string):
		string_number = '0'+string_number
	#end for
	return string_number[0:4]
#end def number_to_4_string


def He_W_xolotl(sim_number=1, number_histories=1000000,
	incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
	IQ0=0, number_layers=100, depth=200.0):

	#Lookup species physical parameters from table
	lookup_table = species_lookup_table()
	Helium = lookup_table.find_species('He', E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	Tungsten = lookup_table.find_species('W', QU=1.0, QUMAX=1.0)

	#Collect species into species_list
	species_list = [Helium, Tungsten]
	number_species = np.size(species_list)

	#Create fractal and print surface file to 1p###.surf
	surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
	surface.print_fractal_surface_file()

	#Define simulation parameters for F-TRIDYN
	simulation_parameters = sim_params(name='He_W', sim_num=sim_number, IFOUT=int(number_histories/20),
		NH=int(number_histories), IDOUT=int(number_histories/5), IQOUT=int(number_histories/5), NCP=number_species, IDREL=1,
		IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
		IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
		DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
		output_pka=0, output_ska=0, output_prj=1)

	#Print simulation parameters to <name><sim_number>.IN
	simulation_parameters.print_input_file()
#end def He_W

def Prj_Tg_xolotl(sim_number=1, number_histories=1000000,
                  incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
                  IQ0=0, number_layers=100, depth=200.0, projectile_name='He',target1_name='W', target2_name='',target3_name='',target4_name=''):

        #Lookup species physical parameters from table       
        lookup_table = species_lookup_table()
        projectile=lookup_table.find_species(projectile_name, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
        target1=lookup_table.find_species(target1_name, QU=1.0, QUMAX=1.0)

        #Collect species into species_list                   
        species_list = [projectile,target1]

        #add target species if other than empty: 
        #currently max of 5 species handled by F-Tridyn 
        additional_targets=[target2_name, target3_name, target4_name]

        for  sp in additional_targets:
                if sp:
                        add_target=lookup_table.find_species(sp, QU=1.0, QUMAX=1.0)
                        species_list.append(add_target)
                

        number_species = np.size(species_list)

        #Create fractal and print surface file to 1p###.surf 
        surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
        surface.print_fractal_surface_file()

        #string of 4-characters required for name
        if (len(projectile_name)+len(target1_name))==4:
                parameterNameString=projectile_name+target1_name #already 4 char, e.g. HeTa
        elif (len(projectile_name)+len(target1_name))==3: 
                parameterNameString=projectile_name+'_'+target1_name #e.g.: He_W  
        else: #len(parameterNameString) ==2
                parameterNameString=projectile_name+'_'+target1_name+'_' #e.g., W_W -> W_W_

        #Define simulation parameters for F-TRIDYN           
        simulation_parameters = sim_params(name=parameterNameString, sim_num=sim_number, IFOUT=int(number_histories/20),
                NH=int(number_histories), IDOUT=int(number_histories/5), IQOUT=int(number_histories/5), NCP=number_species, IDREL=1,
                IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
                IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
                DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
                output_pka=0, output_ska=0, output_prj=1)

        #Print simulation parameters to <name><sim_number>.IN
        simulation_parameters.print_input_file()
#end def Prj_Tg

def beam_and_target(name,beam_species,target_species,sim_number=1,
	number_histories=1E5,incident_energy=100.0,incident_angle=0.0,
	fractal_dimension=1.0,width=200.0,depth=200.0,fluence = 100.0,dataFile='table1.txt'):

	lookup_table = species_lookup_table(dataFile)
	beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	target = lookup_table.find_species(target_species, QU=1.0, QUMAX=1.0)
	species_list = [beam,target]
	number_species = np.size(species_list)

	surface = fractal_surface(FD=fractal_dimension, width=width, iterations=3)
	surface.print_fractal_surface_file()

	simulation_parameters = sim_params(name=name, sim_num=sim_number, IFOUT=number_histories/20,
		NH=number_histories, IDOUT=number_histories, IQOUT=number_histories, NCP=number_species, IDREL=0,
		IQ0=0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=fluence, INEL=1,
		IWC=3, IDIFF=0, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=500,
		DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
		output_pka=0, output_ska=0, output_prj=1)
	simulation_parameters.print_input_file()
#end def beam_and_target

def iterate_SBV(name, beam_species, target_species, incident_energy, incident_angle, depth, SBV_number, SBV_change):
	number_histories = 1E5
	lookup_table = species_lookup_table()
	beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	target = lookup_table.find_species(target_species, QU=1.0, QUMAX=1.0)
	SBV_values = target.SBV * np.arange(1.0-SBV_change,1.0+SBV_change,(2.0*SBV_change / SBV_number))
	print(SBV_values)
	surface = fractal_surface(FD=1.0, width=200, iterations=1)
	surface.print_fractal_surface_file()

	for i in range(SBV_number):
		target.SBV = SBV_values[i]
		species_list = [beam,target]
		number_species = np.size(species_list)
		simulation_parameters = sim_params(name=name, sim_num=i+1, IFOUT=number_histories/20,
			NH=number_histories, IDOUT=number_histories/5, IQOUT=number_histories/5, NCP=number_species, IDREL=0,
			IQ0=0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1.0E-16, INEL=1,
			IWC=3, IDIFF=0, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=500,
			DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
			output_pka=0, output_ska=0, output_prj=1)
		simulation_parameters.print_input_file()

class ftridyn_output_data:
	def __init__(self,name,simulation_number):
		output_file_lines = open(name+simulation_number+'.OUT','r').readlines()

		data_type = []
		species_number = []

		for i in range(len(output_file_lines)):
			data_type.append(output_file_lines[i][0])
			species_number.append(output_file_lines[i][2])
		#end for

		#Projectiles
		self.projectile_species = []
		self.number_projectiles = []
		self.total_projectile_energy = []
		for i in range(data_type.index('P'),data_type.index('P')+data_type.count('P')):
			self.projectile_species.append(int(species_number[i]))
			self.number_projectiles.append(int(output_file_lines[i][4:11].strip()))
			self.total_projectile_energy.append(float(output_file_lines[i][12:22].strip()))
		#end for
		self.average_incident_energy = np.sum(self.total_projectile_energy) / np.sum(self.number_projectiles)

		#Backscattering
		self.backscattered_species = []
		self.number_backscattered = []
		self.total_backscattered_energy = []
		for i in range(data_type.index('B'),data_type.index('B')+data_type.count('B')):
			self.backscattered_species.append(int(species_number[i]))
			self.number_backscattered.append(int(output_file_lines[i][4:11].strip()))
			self.total_backscattered_energy.append(float(output_file_lines[i][12:22].strip()))
		#end for

		#Transmission
		self.transmitted_species = []
		self.number_transmitted = []
		self.total_transmitted_energy = []
		for i in range(data_type.index('T'),data_type.index('T')+data_type.count('T')):
			self.transmitted_species.append(int(species_number[i]))
			self.number_transmitted.append(int(output_file_lines[i][4:11].strip()))
			self.total_transmitted_energy.append(float(output_file_lines[i][12:22].strip()))
		#end for

		#Sputtering
		self.sputtered_species = []
		self.number_sputtered = []
		self.total_sputtered_energy = []
		for i in range(data_type.index('S'),data_type.index('S')+data_type.count('S')):
			self.sputtered_species.append(int(species_number[i]))
			self.number_sputtered.append(int(output_file_lines[i][4:11].strip()))
			self.total_sputtered_energy.append(float(output_file_lines[i][12:22].strip()))
		#end for

		#Frenkel Pairs
		self.frenkel_pair_species = []
		self.number_frenkel_pairs = []
		for i in range(data_type.index('F'),data_type.index('F')+data_type.count('F')):
			self.frenkel_pair_species.append(int(species_number[i]))
			self.number_frenkel_pairs.append(float(output_file_lines[i][4:13].strip()))
		#end for
	#end def __init__

	def calculate_total_sputtering_yield(self):
		number_sputtered = np.sum(self.number_sputtered)
		sputtering_yield = float(number_sputtered) / float(np.sum(self.number_projectiles))
		return sputtering_yield
	#end def calculate_total_sputtering_yield
	def calculate_total_reflection_yield(self):
		number_reflected = np.sum(self.number_backscattered)
		reflection_yield = float(number_reflected) / float(np.sum(self.number_projectiles))
		return reflection_yield

def number_to_4_string(number):
	string_number = str(number)
	#todo: try except blocks to capture error
	length_string = len(string_number)
	for i in range(4-length_string):
		string_number = '0'+string_number
	#end for
	return string_number[0:4]
#end def number_to_4_string

def plot_sputtered_angular_distributions(simulation_name,nAbins,plotsOn = 0):
    splst = np.loadtxt(simulation_name+'SPLST.DAT')
    print('splist size ',splst.size)
    #n = np.zeros(shape=(50))
    #bins = np.zeros(shape=(51))
    if splst.size > 10 :
        cos_angles = splst[:,6:9]
        alpha0 = cos_angles[:,0]
        beta0 = cos_angles[:,1]
        gamma0 = cos_angles[:,2]
        angles = np.arccos(cos_angles)*180.0 / np.pi
        theta = np.arctan2(gamma0,beta0)*180.0/np.pi
        alpha = angles[:,0]
        beta = angles[:,1]
        gamma = angles[:,2]
        nX,binsX,patches = plt.hist(alpha,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.close()    
            plt.figure(1)
            plt.title(simulation_name+' Alpha Angle Distribution')
            plt.xlabel('alpha [deg]')
            plt.savefig(simulation_name+'alpha_dist.png')
        
        nY,binsY,patches = plt.hist(theta,bins=nAbins,range=(-180.0,180.0))
        if plotsOn :
            plt.figure(2)
            plt.title(simulation_name+' Beta Angle Distribution')
            plt.xlabel('beta [deg]')
            plt.savefig(simulation_name+'beta_dist.png')
        
        #plt.figure(3)
        #nZ,binsZ,patches = plt.hist(gamma,bins=nAbins,range=(0.0,90.0))
        #if plotsOn :
        #    plt.title(simulation_name+' Gamma Angle Distribution')
        #    plt.xlabel('gamma [deg]')
        #    plt.savefig(simulation_name+'gamma_dist.png')
    else :
        nX = np.zeros(shape=(nAbins))
        binsX = np.zeros(shape=(nAbins+1))
        nY = nX
        #nZ = nX
        binsY = binsX
        #binsZ = binsX

    return (splst.size, nX, binsX,nY, binsY)

def plot_sputtered_energy_distributions(simulation_name,nEbins,maxE=100.0, plotsOn=0):
    splst = np.loadtxt(simulation_name+'SPLST.DAT')
    print('splist size ',splst.size)
    if splst.size > 10 :
        en = splst[:,2]
        plt.close()    
        plt.figure(1)
        n,bins,patches = plt.hist(en,bins=nEbins,range=(0.0,maxE))
        if plotsOn :
            plt.title(simulation_name+' Energy Distribution')
            plt.xlabel('Energy [eV]')
            plt.savefig(simulation_name+'alpha_dist.png')
    else :
        n = np.zeros(shape=(nEbins))
        bins = np.zeros(shape=(nEbins+1))

    return (splst.size, n, bins)

def plot_reflected_angular_distributions(simulation_name,nAbins,plotsOn = 0):
    splst = np.loadtxt(simulation_name+'RFLST.DAT')
    print('rflist size ',splst.size)
    #n = np.zeros(shape=(50))
    #bins = np.zeros(shape=(51))
    if splst.size > 10 :
        cos_angles = splst[:,6:9]
        alpha0 = cos_angles[:,0]
        beta0 = cos_angles[:,1]
        gamma0 = cos_angles[:,2]
        angles = np.arccos(cos_angles)*180.0 / np.pi
        theta = np.arctan2(gamma0,beta0)*180.0/np.pi

        alpha = angles[:,0]
        beta = angles[:,1]
        gamma = angles[:,2]
        nX,binsX,patches = plt.hist(alpha,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.close()    
            plt.figure(1)
            plt.title(simulation_name+' Alpha Angle Distribution')
            plt.xlabel('alpha [deg]')
            plt.savefig(simulation_name+'alpha_dist.png')
        
        nY,binsY,patches = plt.hist(theta,bins=nAbins,range=(-180.0,180.0))
        if plotsOn :
            plt.figure(2)
            plt.title(simulation_name+' Beta Angle Distribution')
            plt.xlabel('beta [deg]')
            plt.savefig(simulation_name+'beta_dist.png')
        
        #plt.figure(3)
        #nZ,binsZ,patches = plt.hist(gamma,bins=nAbins,range=(0.0,90.0))
        #if plotsOn :
        #    plt.title(simulation_name+' Gamma Angle Distribution')
        #    plt.xlabel('gamma [deg]')
        #    plt.savefig(simulation_name+'gamma_dist.png')
    else :
        nX = np.zeros(shape=(nAbins))
        binsX = np.zeros(shape=(nAbins+1))
        nY = nX
        #nZ = nX
        binsY = binsX
        #binsZ = binsX

    return (splst.size, nX, binsX,nY, binsY)

def plot_reflected_energy_distributions(simulation_name,nEbins=500,maxE=1000.0, plotsOn=0):
    splst = np.loadtxt(simulation_name+'RFLST.DAT')
    print('rflist size ',splst.size)
    if splst.size > 10 :
        en = splst[:,2]
        plt.close()    
        plt.figure(1)
        print('en bins range',en, nEbins, maxE)
        n,bins,patches = plt.hist(en,bins=nEbins,range=(0.0,maxE))
        if plotsOn :
            plt.title(simulation_name+' Energy Distribution')
            plt.xlabel('Energy [eV]')
            plt.savefig(simulation_name+'alpha_dist.png')
    else :
        n = np.zeros(shape=(nEbins))
        bins = np.zeros(shape=(nEbins+1))

    return (splst.size, n, bins)

def plot_incident_energy_distribution(simulation_name, simulation_number):
    ed1 = np.loadtxt(simulation_name+number_to_4_string(simulation_number)+'.ED1')
    plt.close()
    plt.figure(1)
    plt.plot(ed1[:,0],ed1[:,1])
    plt.title(simulation_name+number_to_4_string(simulation_number)+' Incident Energy Distribution')
    plt.xlabel('E [eV]')
    plt.savefig(simulation_name+number_to_4_string(simulation_number)+'incident_E_dist.png')

def plot_implantation_profile(simulation_name,plotsOn=0):
    prj = np.loadtxt(simulation_name+'DUMPPRJ.dat')
    x = prj[:,2]
    if plotsOn :
        plt.close()
        plt.figure(1)
        plt.hist(x,bins=100)
        plt.title(simulation_name+' Implantation Profile')
        plt.xlabel('x [Ang]')
        plt.savefig(simulation_name+'implantation.png')

def plot_sputtering_yield(simulation_name,simulation_numbers):
	simulation_output = []
	energy = []
	sputtering_yield = []
	for i in range(np.size(simulation_numbers)):
		simulation_output.append(ftridyn_output_data(simulation_name, number_to_4_string(simulation_numbers[i])))
		energy.append(simulation_output[i].average_incident_energy)
		sputtering_yield.append(simulation_output[i].calculate_total_sputtering_yield())
		print(energy[i],sputtering_yield[i])
	plt.loglog(energy,sputtering_yield)
	plt.title(simulation_name)
	plt.show()
def func1(path,E,a,r,d,specNum):
    print(path)
    cwd=os.getcwd()
    os.mkdir(path)
    os.chdir(path)
    #cwd1=os.getcwd()
    #print('cwd1 ',cwd1)
    beam = d['beam'][specNum]
    name1 = ''
    name2 = ''
    if(len(beam)>1):
        name1 = beam
    else:
        name1 = beam+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']

    beam_and_target(name1+name2,beam,d['target'],sim_number=1,number_histories=d['nH'], incident_energy=E,depth=200.0,incident_angle=a,fluence=1.0E-16,dataFile=d['data'])
    #p = subprocess.Popen([d['exe'],name1+name2+'0001.IN'],cwd=cwd+'/'+path,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #stdoutdata, stderrdata = p.communicate()
    #returnCode = p.wait()
    #try:
    #    fileLog = open(path+'/log.txt','w') 
    #    fileErr = open(path+'/logErr.txt','w')
    #    fileIn = open(path+'/'+name1+name2+'0001.IN')
    #    p = subprocess.check_call([d['exe']],stdin=fileIn,cwd=cwd+'/'+path,stdout=fileLog, stderr=fileErr)
    #    fileLog.close()
    #    fileErr.close()
    #except subprocess.CalledProcessError as e:
    #    sys.exit("'ls' failed, returned code %d (check 'errors.txt')" \
    #            % (e.returncode))
    #ftridyn.ftridyn_cpmi(name1+name2+'0001.IN')
    #returnCode=0
    if d['use_exe']:
        try:
            #original = sys.stdout
            f1 =  open('log.txt', 'w')
            f1.close()
            f2 =  open('logErr.txt', 'w')
            f2.close()
            null_fds = [os.open('log.txt',os.O_RDWR), os.open('logErr.txt',os.O_RDWR)]
            # save the current file descriptors to a tuple
            save = os.dup(1), os.dup(2)
            # put /dev/null fds on 1 and 2
            os.dup2(null_fds[0], 1)
            os.dup2(null_fds[1], 2)
            #ftridyn.ftridyn_cpmi(name1+name2+'0001.IN')
            ff = subprocess.check_output(d['exe']+' '+name1+name2+'0001.IN', shell=True)
            #ff = subprocess.check_output('/Users/tyounkin/Code/ftridyn2/src/shell_FtridynGITR.sh '+ name1+name2+'0001.IN', shell=True)
            # restore file descriptors so I can print the results
            os.dup2(save[0], 1)
            os.dup2(save[1], 2)
            # close the temporary fds
            os.close(null_fds[0])
            os.close(null_fds[1])
            os.close(save[0])
            os.close(save[1])
            print('ran tridyn')
            returnCode=0
        #sys.stdout = original
        #print 'ftridyn print', mystdout.getValue()
        except subprocess.CalledProcessError as e:
            sys.stderr.write("'ls' failed, returned code %d (check 'errors.txt')\n" % (e.returncode))
#        except:
#            print("FTRIDYN ERROR!")
            returnCode=0
    else:    
        try:
            #original = sys.stdout
            f1 =  open('log.txt', 'w')
            f1.close()
            f2 =  open('logErr.txt', 'w')
            f2.close()
            null_fds = [os.open('log.txt',os.O_RDWR), os.open('logErr.txt',os.O_RDWR)]
            # save the current file descriptors to a tuple
            save = os.dup(1), os.dup(2)
            # put /dev/null fds on 1 and 2
            os.dup2(null_fds[0], 1)
            os.dup2(null_fds[1], 2)
            #ftridyn.ftridyn_cpmi(name1+name2+'0001.IN')
            # restore file descriptors so I can print the results
            os.dup2(save[0], 1)
            os.dup2(save[1], 2)
            # close the temporary fds
            os.close(null_fds[0])
            os.close(null_fds[1])
            os.close(save[0])
            os.close(save[1])
            print('ran tridyn')
            returnCode=0
        #sys.stdout = original
        #print 'ftridyn print', mystdout.getValue()
        except:
            print("FTRIDYN ERROR!")
            returnCode=0
    #print('path, returncode', path, p)
    #file = open(path+'/log.txt','w') 
    #file.write(str(stdoutdata)) 
    #file.close() 
    #file = open(path+'/logErr.txt','w') 
    #file.write(str(stderrdata)) 
    #file.close() 
    #ftridyn.ftridyn_cpmi('He_W0001.IN')
    ##time.sleep(60)
    #print('finished simulation')
    #sim = analyze_ftridyn_simulations.ftridyn_output_data('He_W','0001')
    #Y = sim.calculate_total_sputtering_yield()
    #print('Energy Yield ',E, Y)
    os.chdir(cwd)
    #cwd=os.getcwd()
    analyzeFTrun(path,d,specNum)
    print('cwd ',cwd)
    return returnCode
def analyzeFTrun(pathString,d,specNum):
    nEgrid = d['nEdist']
    maxE = d['maxEdist']
    nEgrid_ref = d['nEdist_ref']
    maxE_ref = d['maxEdist_ref']
    nAgrid = d['nAdist']
    beam = d['beam'][specNum]
    name1 = ''
    name2 = ''
    if(len(beam)>1):
        name1 = beam
    else:
        name1 = beam+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']
    ffilename = name1+name2
    spyl_file = ffilename +'SPYL.DAT'
    print('doing analysis on ' , pathString+"/"+spyl_file)
    WW=ftridyn_output_data(pathString+"/"+ffilename,'0001')
    thisSpyl = WW.calculate_total_sputtering_yield()
    print('sputtering yield', thisSpyl)
    thisRefyl = WW.calculate_total_reflection_yield()
    np.savetxt(pathString+"/"+"YR.out",[thisSpyl,thisRefyl])
    print('reflection yield', thisRefyl)
    nP, nX, binsX,nY,binsY = plot_sputtered_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nX.out",nX)
    np.savetxt(pathString+"/"+"nY.out",nY)
    #np.savetxt(pathString+"/"+"nZ.out",nZ)
    nP, nX, binsX,nY,binsY = plot_reflected_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nXref.out",nX)
    np.savetxt(pathString+"/"+"nYref.out",nY)
    #np.savetxt(pathString+"/"+"nZref.out",nZ)
    nPenergy, nenergy, binsenergy = plot_sputtered_energy_distributions(pathString+"/"+ffilename,nEgrid,maxE)
    np.savetxt(pathString+"/"+"energy.out",nenergy)
    nPenergyRef, nenergyRef, binsenergyRef = plot_reflected_energy_distributions(pathString+"/"+ffilename,nEgrid_ref,maxE_ref)
    np.savetxt(pathString+"/"+"energyRef.out",nenergyRef)

def func2(path,E,a,r):
     print('Path ear ', E, a , r)
class Work():
    def __init__(self,work_items):
        self.work_items = work_items[:] 
 
    def get_next_item(self):
        if len(self.work_items) == 0:
            return None
        return self.work_items.pop()
def master(nList):
    indxList = list(range(1,nList+1))
    all_data = []
    np = MPI.COMM_WORLD.Get_size()
    current_work = Work(indxList) 
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    size = min(np,nList)
    for i in range(1, size): 
        anext = current_work.get_next_item() 
        if not anext: break
        comm.send(anext-1, dest=i, tag=WORKTAG)
    
    while 1:
        anext = current_work.get_next_item()
        if not anext: break
        data = comm.recv(None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        all_data.append(data)
        comm.send(anext-1, dest=status.Get_source(), tag=WORKTAG)
 
    for i in range(1,size):
        data = comm.recv(None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
        all_data.append(data)
    
    for i in range(1,np):
        comm.send(None, dest=i, tag=DIETAG)
     
    return all_data
        
    
def slave(func,pathList,eList,aList,rList,d):
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(None, source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag(): break
        else:
            print('rank %d (total %d) received data %i' % ( rank, size,int(data)) )
        specNum = int(data)%d['nS']
        complete = func(pathList[data],eList[data],aList[data],rList[data],d,specNum)
        comm.send(complete, dest=0)

def main(func,argv):
    print('Argument List:', str(argv))
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    comm = MPI.COMM_WORLD
    exit=0
    print('I am  %s rank %d (total %d)' % (name, rank, size) )
    loadDict = False
    dictName = None
    if(rank ==0):
        start_time = time.time()
        try:
            opts, args = getopt.getopt(argv, "d", ["dictionary="])
        except getopt.GetoptError as err:
            # print help information and exit:
            print(str(err))  # will print something like "option -a not recognized"
            #usage()
            exit=1
            #sys.exit(2)
        output = None
        verbose = False
        print("opts ", opts)
        if(exit==0):
            for o, a in opts:
                if o in ("-d", "--dictionary"):
                    #usage()
                    #sys.exit()
                    print("dictionary " , a)
                    loadDict = True
                    dictName = a
                else:
                    assert False, "unhandled option"
                    exit=1

    exit = comm.bcast(exit,root=0)
    if(exit):
        MPI.COMM_WORLD.Abort()
        sys.exit()
    pathList = []
    eList = []
    aList = []
    rList = []
    #d = {}
    d = defaultdict(list)
    
    if rank == 0: # Master
        print('master rank',rank)
        if loadDict:
            f = open(dictName, 'r')   # 'r' for reading; can be omitted
            d = pickle.load(f)         # load file content as mydict
            f.close()
        else:
            with open("ftMPI.in") as f:
                for line in f:
                    split = line.split()
                    if(split[0] == 'beam' or split[0]=='target' or split[0]=='Escale' or split[0]=='exe' or split[0]=='data'):
                        for j in range(1,len(split)):
                            d[split[0]].append(split[j])
                    elif(split[0] == 'nE' or split[0] == 'nA' or split[0]=='nR' or split[0]=='nEdist' or split[0]=='nEdist_ref' or split[0]=='nAdist' or split[0]=='nH' or split[0]=='use_exe'):
                        d[split[0]] = int(split[1])
                    else:
                        d[split[0]] = float(split[1])
        print(d)
        d['Escale'] = d['Escale'][0]
        d['target'] = d['target'][0]
        d['exe'] = d['exe'][0]
        d['data'] = d['data'][0]
        if(d['Escale'] == 'log'):
            energy = np.logspace(d['energy_start'],d['energy_end'],d['nE'])
        else:
            energy = np.linspace(d['energy_start'],d['energy_end'],d['nE'])
        angle = np.linspace(d['angle_start'],d['angle_end'],d['nA'])
        roughness = np.linspace(d['roughness_start'],d['roughness_end'],d['nR'])
        d['nS'] = len(d['beam'])
        print('nS', d['nS'], d['beam'])
        for i in range(d['nE']):
            for j in range(d['nA']):
                for k in range(d['nR']):
                    for ii in range(d['nS']):
                        pathString = "FTRIDYN_"+d['beam'][ii]+"_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                        pathList.append(pathString)
                        eList.append(energy[i])
                        aList.append(angle[j])
                        rList.append(roughness[k])
        print(pathList)
        f = open('pathList.pkl', 'wb')
        pickle.dump(pathList,f)
        f.close() 
    comm.Barrier()
    pathList = comm.bcast(pathList,root=0)
    eList = comm.bcast(eList,root=0)
    aList = comm.bcast(aList,root=0)
    rList = comm.bcast(rList,root=0)
    d = comm.bcast(d,root=0)
    comm.Barrier()
    all_data = []
    if rank == 0: # Master
        print('master',rank)
        all_data = master(len(pathList))
    else: # Any slave
        print('slave rank',rank)
        slave(func,pathList,eList,aList,rList,d)
    
    print('Task waiting at barrier (rank %d)' % (rank) )
    comm.Barrier()
    print('Task completed (rank %d)' % (rank) )
    if(rank ==0):
        print('returns ' , all_data)
    
    if rank == 0:
        nEgrid = d['nEdist']
        maxE = d['maxEdist']
        nEgrid_ref = d['nEdist_ref']
        maxE_ref = d['maxEdist_ref']
        nAgrid = d['nAdist']
        sputt = np.zeros(shape=(d['nS'],len(energy),len(angle)))
        refl = np.zeros(shape=(d['nS'],len(energy),len(angle)))
        eDistEgrid = np.linspace(0.0,maxE-maxE/nEgrid,nEgrid) 
        eDistEgridRef = np.linspace(0.0,maxE_ref-maxE_ref/nEgrid_ref,nEgrid_ref) 
        phiGrid = np.linspace(0.0,90.0-90.0/nAgrid,nAgrid) 
        thetaGrid = np.linspace(0.0,180.0-180.0/nAgrid,nAgrid) 
        cosXDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosYDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        #cosZDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosXDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosYDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        #cosZDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        eDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nEgrid)) 
        eDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nEgrid_ref))
        sputtSelf = np.zeros(shape=(len(energy),len(angle)))
        reflSelf = np.zeros(shape=(len(energy),len(angle)))
        cosXDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        #cosZDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosXDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        #cosZDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        eDistributionSelf = np.zeros(shape=(len(energy),len(angle),nEgrid)) 
        eDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nEgrid_ref))
        totalIndex=0
        for i in range(len(energy)):
            for j in range(len(angle)):
                for ii in range(d['nS']):
                    specNum = ii%d['nS']
                    pathString = pathList[totalIndex]
                    yr = np.loadtxt(pathString+"/"+"YR.out", dtype='float')
                    sputt[specNum,i,j] = yr[0]
                    refl[specNum,i,j] = yr[1]
                    nX = np.loadtxt(pathString+"/"+"nX.out", dtype='float')
                    cosXDistribution[specNum,i,j,:] = nX
                    nY = np.loadtxt(pathString+"/"+"nY.out", dtype='float')
                    cosYDistribution[specNum,i,j,:] = nY
                    #nZ = np.loadtxt(pathString+"/"+"nZ.out", dtype='float')
                    #cosZDistribution[specNum,i,j,:] = nZ
                    nXref = np.loadtxt(pathString+"/"+"nXref.out", dtype='float')
                    cosXDistributionRef[specNum,i,j,:] = nXref
                    nYref = np.loadtxt(pathString+"/"+"nYref.out", dtype='float')
                    cosYDistributionRef[specNum,i,j,:] = nYref
                    #nZref = np.loadtxt(pathString+"/"+"nZref.out", dtype='float')
                    #cosZDistributionRef[specNum,i,j,:] = nZref
                    nenergy = np.loadtxt(pathString+"/"+"energy.out", dtype='float')
                    eDistribution[specNum,i,j,:] = nenergy
                    nenergyRef = np.loadtxt(pathString+"/"+"energyRef.out", dtype='float')
                    eDistributionRef[specNum,i,j,:] = nenergyRef
                    totalIndex = totalIndex+1
	#for i in range(d['nS']):
        if d['beam'][-1] == d['target']:
            nS_background = d['nS']-1
            save_self = 1
        else:
            nS_background = d['nS']
            save_self = 0

        rootgrp = netCDF4.Dataset("ftridynBackground"+".nc", "w", format="NETCDF4")
        ne = rootgrp.createDimension("nE", len(energy))
        na = rootgrp.createDimension("nA", len(angle))
        ns = rootgrp.createDimension("nS", nS_background)
        nedistgrid = rootgrp.createDimension("nEdistBins", nEgrid)
        nedistgridref = rootgrp.createDimension("nEdistBinsRef", nEgrid_ref)
        nadistgrid = rootgrp.createDimension("nAdistBins", nAgrid)
        spyld = rootgrp.createVariable("spyld","f8",("nS","nE","nA"))
        rfyld = rootgrp.createVariable("rfyld","f8",("nS","nE","nA"))
        ee = rootgrp.createVariable("E","f8",("nE"))
        aa = rootgrp.createVariable("A","f8",("nA"))
        cosxdist = rootgrp.createVariable("cosXDist","f8",("nS","nE","nA","nAdistBins"))
        cosydist = rootgrp.createVariable("cosYDist","f8",("nS","nE","nA","nAdistBins"))
        #coszdist = rootgrp.createVariable("cosZDist","f8",("nS","nE","nA","nAdistBins"))
        cosxdistref = rootgrp.createVariable("cosXDistRef","f8",("nS","nE","nA","nAdistBins"))
        cosydistref = rootgrp.createVariable("cosYDistRef","f8",("nS","nE","nA","nAdistBins"))
        #coszdistref = rootgrp.createVariable("cosZDistRef","f8",("nS","nE","nA","nAdistBins"))
        edist = rootgrp.createVariable("energyDist","f8",("nS","nE","nA","nEdistBins"))
        edistref = rootgrp.createVariable("energyDistRef","f8",("nS","nE","nA","nEdistBinsRef"))
        edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins")) 
        edistegridref = rootgrp.createVariable("eDistEgridRef","f8",("nEdistBinsRef")) 
        phigrid = rootgrp.createVariable("phiGrid","f8",("nAdistBins")) 
        thetagrid = rootgrp.createVariable("thetaGrid","f8",("nAdistBins")) 
        ee[:] = energy
        aa[:] = angle
        edistegrid[:] = eDistEgrid
        edistegridref[:] = eDistEgridRef
        phigrid[:] = phiGrid
        thetagrid[:] = thetaGrid
        spyld[:] = sputt[0:(nS_background),:,:]
        rfyld[:] = refl[0:(nS_background),:,:]
        cosxdist[:] = cosXDistribution[0:(nS_background),:,:,:]
        cosydist[:] = cosYDistribution[0:(nS_background),:,:,:]
        #coszdist[:] = cosZDistribution[0:d['nS'],:,:,:]
        cosxdistref[:] = cosXDistributionRef[0:(nS_background),:,:,:]
        cosydistref[:] = cosYDistributionRef[0:(nS_background),:,:,:]
        #coszdistref[:] = cosZDistributionRef[0:d['nS'],:,:,:]
        edist[:] = eDistribution[0:(nS_background),:,:,:]
        edistref[:] = eDistributionRef[0:(nS_background),:,:,:]
        rootgrp.close()
        
        if save_self:
            rootgrp = netCDF4.Dataset("ftridynSelf"+".nc", "w", format="NETCDF4")
            ne = rootgrp.createDimension("nE", len(energy))
            na = rootgrp.createDimension("nA", len(angle))
            nedistgrid = rootgrp.createDimension("nEdistBins", nEgrid)
            nedistgridref = rootgrp.createDimension("nEdistBinsRef", nEgrid_ref)
            nadistgrid = rootgrp.createDimension("nAdistBins", nAgrid)
            spyld = rootgrp.createVariable("spyld","f8",("nE","nA"))
            rfyld = rootgrp.createVariable("rfyld","f8",("nE","nA"))
            ee = rootgrp.createVariable("E","f8",("nE"))
            aa = rootgrp.createVariable("A","f8",("nA"))
            cosxdist = rootgrp.createVariable("cosXDist","f8",("nE","nA","nAdistBins"))
            cosydist = rootgrp.createVariable("cosYDist","f8",("nE","nA","nAdistBins"))
            #coszdist = rootgrp.createVariable("cosZDist","f8",("nE","nA","nAdistBins"))
            cosxdistref = rootgrp.createVariable("cosXDistRef","f8",("nE","nA","nAdistBins"))
            cosydistref = rootgrp.createVariable("cosYDistRef","f8",("nE","nA","nAdistBins"))
            #coszdistref = rootgrp.createVariable("cosZDistRef","f8",("nE","nA","nAdistBins"))
            edist = rootgrp.createVariable("energyDist","f8",("nE","nA","nEdistBins"))
            edistref = rootgrp.createVariable("energyDistRef","f8",("nE","nA","nEdistBinsRef"))
            edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins")) 
            edistegridref = rootgrp.createVariable("eDistEgridRef","f8",("nEdistBinsRef")) 
            phigrid = rootgrp.createVariable("phiGrid","f8",("nAdistBins")) 
            thetagrid = rootgrp.createVariable("thetaGrid","f8",("nAdistBins")) 
            phigrid[:] = phiGrid
            thetagrid[:] = thetaGrid
            ee[:] = energy
            aa[:] = angle
            edistegrid[:] = eDistEgrid
            edistegridref[:] = eDistEgridRef
            spyld[:] = sputt[-1,:,:]
            rfyld[:] = refl[-1,:,:]
            cosxdist[:] = cosXDistribution[-1,:,:,:]
            cosydist[:] = cosYDistribution[-1,:,:,:]
            #coszdist[:] = cosZDistribution[-1,:,:,:]
            cosxdistref[:] = cosXDistributionRef[-1,:,:,:]
            cosydistref[:] = cosYDistributionRef[-1,:,:,:]
            #coszdistref[:] = cosZDistributionRef[-1,:,:,:]
            edist[:] = eDistribution[-1,:,:,:]
            edistref[:] = eDistributionRef[-1,:,:,:]
            rootgrp.close()

        exec_time = time.time()
        print("Execution of FTRIDYN Cases took --- %s seconds ---" % (exec_time - start_time))

if __name__ == "__main__":
    print('Argument List:', str(sys.argv[1:]))
    #print('Using ftridyn library from ', ftridyn.__file__)
    main(func1,sys.argv[1:])
    #main(func2)
