from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
        angles = np.arccos(cos_angles)*180.0 / np.pi

        alpha = angles[:,0]
        beta = angles[:,1]
        gamma = angles[:,2]
        plt.close()    
        plt.figure(1)
        nX,binsX,patches = plt.hist(alpha,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Alpha Angle Distribution')
            plt.xlabel('alpha [deg]')
            plt.savefig(simulation_name+'alpha_dist.png')
        
        plt.figure(2)
        nY,binsY,patches = plt.hist(beta,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Beta Angle Distribution')
            plt.xlabel('beta [deg]')
            plt.savefig(simulation_name+'beta_dist.png')
        
        plt.figure(3)
        nZ,binsZ,patches = plt.hist(gamma,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Gamma Angle Distribution')
            plt.xlabel('gamma [deg]')
            plt.savefig(simulation_name+'gamma_dist.png')
    else :
        nX = np.zeros(shape=(nAbins))
        binsX = np.zeros(shape=(nAbins+1))
        nY = nX
        nZ = nX
        binsY = binsX
        binsZ = binsX

    return (splst.size, nX, binsX,nY, binsY,nZ, binsZ)

def plot_sputtered_energy_distributions(simulation_name,nEbins, plotsOn=0):
    splst = np.loadtxt(simulation_name+'SPLST.DAT')
    print('splist size ',splst.size)
    if splst.size > 10 :
        en = splst[:,2]
        plt.close()    
        plt.figure(1)
        n,bins,patches = plt.hist(en,bins=nEbins,range=(0.0,100.0))
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
        angles = np.arccos(cos_angles)*180.0 / np.pi

        alpha = angles[:,0]
        beta = angles[:,1]
        gamma = angles[:,2]
        plt.close()    
        plt.figure(1)
        nX,binsX,patches = plt.hist(alpha,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Alpha Angle Distribution')
            plt.xlabel('alpha [deg]')
            plt.savefig(simulation_name+'alpha_dist.png')
        
        plt.figure(2)
        nY,binsY,patches = plt.hist(beta,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Beta Angle Distribution')
            plt.xlabel('beta [deg]')
            plt.savefig(simulation_name+'beta_dist.png')
        
        plt.figure(3)
        nZ,binsZ,patches = plt.hist(gamma,bins=nAbins,range=(0.0,90.0))
        if plotsOn :
            plt.title(simulation_name+' Gamma Angle Distribution')
            plt.xlabel('gamma [deg]')
            plt.savefig(simulation_name+'gamma_dist.png')
    else :
        nX = np.zeros(shape=(nAbins))
        binsX = np.zeros(shape=(nAbins+1))
        nY = nX
        nZ = nX
        binsY = binsX
        binsZ = binsX

    return (splst.size, nX, binsX,nY, binsY,nZ, binsZ)

def plot_reflected_energy_distributions(simulation_name,nEbins=50, plotsOn=0):
    splst = np.loadtxt(simulation_name+'RFLST.DAT')
    print('rflist size ',splst.size)
    if splst.size > 10 :
        en = splst[:,2]
        plt.close()    
        plt.figure(1)
        n,bins,patches = plt.hist(en,bins=nEbins,range=(0.0,100.0))
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

if __name__ == '__main__':
	plot_sputtering_yield('BeBe',np.arange(1,65))
