#python library tools for gitr

import sys,os
import numpy as np
import io,libconf

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# imports within other functions
# from netCDF4 import Dataset

_VERSION = 1.0
_LAST_UPDATE = 'Aug. 10. 2018'

# ----------------------------------------------------------------------------------------
# general figure font configuration, global change throughout python session
FONTSIZE = 18
if 'Anaconda' in sys.version:
	plt.rcParams['font.sans-serif'] = 'Arial'	# anaconda
	plt.rcParams['font.serif'] = 'Times'	# anaconda
else:
	plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Nimbus Sans L', 'Liberation Sans', 'DejaVu Sans']
	plt.rcParams['font.serif'] = ['Times', 'Times New Roman', 'Nimbus Roman No9 L', 'DejaVu Serif']
FONT = {'family' : 'sans-serif',
		'weight' : 'normal', # normal, bold
		'size'   : FONTSIZE} 
plt.rc('font', **FONT)

# ----------------------------------------------------------------------------------------
# --- gitr Class -------------------------------------------------------------------------

class gitr:
	"""
	Load, Post-process and Plot GITR output
	"""

	def __init__(self, path = './', impurity = 'W', show = False):
		"""
		path (string)      Path to top folder with input and output folders, default './'
		impurity (string)  Main impurity: 'W' (default), ???
		show (bool)        True: make all available plots
		"""
		if not (path[-1] == '/'): path += '/'
		self.path = path
		
		self.config = {}
		# load input configurations
		self.InputFiles = ['gitrInput','gitrGeometry','particleSource','answer']
		for key in self.InputFiles:
			file = self.path + 'input/' + key + '.cfg'
			if os.path.isfile(file): 
				with io.open(file) as f: self.config[key] = libconf.load(f)
		
		self.data = {}
		# load NetCDF output files
		self.files = ['history','positions','spec','surface','particleSource','forces']
		for key in self.files:
			file = self.path + 'output/' + key + '.nc'
			if os.path.isfile(file): self.data[key] = openNCfile(file)
		
		# set impurity
		self.impurity = impurity
		if self.impurity in ['W','w','tungsten','Tungsten']:
			self.Z = 74
			self.M = 184
		
		self.Pmass = 1.66e-27
		self.eCharge = 1.602e-19
		
		# filter NaN
		if self.data.has_key('positions'): 
			self.NotNan = (-np.isnan(self.data['positions']['x'])) & (-np.isnan(self.data['positions']['y'])) & (-np.isnan(self.data['positions']['z']))
		
		if show: self.plotAll()
		
	
	def plotAll(self, N = None):
		"""
		Make all plots 
		N (int)  number of orbits to plot in 3D, default: all
		"""
		self.plotParticles()
		self.plotParticleHistograms()
		self.plotOrbits(N)
		self.plotDensity()
		self.plotEnergyDist()
		
	
	def plotParticles(self):
		"""
		Plot final particle positions
		"""
		if not self.data.has_key('positions'): 
			print 'positions.nc not available'
			return
			
		x = self.data['positions']['x']
		y = self.data['positions']['y']
		z = self.data['positions']['z']
		hit = self.data['positions']['hitWall'] > 0
		notHit = self.data['positions']['hitWall'] == 0
		r = np.sqrt(x**2 + y**2)
				
		fig = plt.figure(figsize = (9,6))
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(x[notHit],y[notHit],z[notHit], c = 'b')
		ax.scatter(x[hit],y[hit],z[hit], c = 'r')
		plt.title('Final ' + self.impurity + ' Particle Positions', fontsize = FONTSIZE)
		ax.view_init(azim = 225)
		ax.ticklabel_format(style = 'sci',useOffset = True, useMathText = True, scilimits = (2,2))
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')

		plt.figure(figsize = (9,6))
		plt.scatter(r[notHit],z[notHit], c = 'b')
		plt.scatter(r[hit],z[hit], c = 'r')
		plt.title('Final ' + self.impurity + ' Particle Positions', fontsize = FONTSIZE)
		plt.xlabel('r [m]')
		plt.ylabel('z [m]')


	def plotParticleHistograms(self):
		"""
		Plot particle, energy and charge deposition patterns on the target surface
		"""
		if not self.data.has_key('positions'): 
			print 'positions.nc not available'
			return

		x = self.data['positions']['x']
		y = self.data['positions']['y']
		z = self.data['positions']['z']
		hit = (self.data['positions']['hitWall'] > 0) & self.NotNan
		
		fig = plt.figure(figsize = (10,6))
		Nx,Ny = 30,30
 		H,yedges,xedges = np.histogram2d(y[hit], x[hit], bins = [Ny,Nx])   # vertical, horizontal
 		X,Y = np.meshgrid(xedges, yedges)
 		cs = plt.imshow(np.log10(H), extent = [X.min(), X.max(), Y.min(), Y.max()], cmap = 'plasma', origin = 'lower', aspect = 'auto', interpolation = 'nearest')
		#cs = plt.pcolormesh(X, Y, H, cmap = 'plasma')
		plt.xlabel('x [m]')
		plt.ylabel('y [m]')
		C = plt.colorbar(cs, pad = 0.01, format = '%.3g')
		C.set_label('Number of Particles', rotation = 270, size = FONTSIZE, va = 'bottom')
		C.ax.set_yticklabels(['10$^{' + item.get_text() + '}$' for item in C.ax.get_yticklabels()])
		plt.title('Target Erosion: Particles', fontsize = FONTSIZE)


		xs = x[hit]
		ys = y[hit]
		energies = np.zeros((Ny,Nx))
		charges = np.zeros((Ny,Nx))

		vx = self.data['positions']['vx']
		vy = self.data['positions']['vy']
		vz = self.data['positions']['vz']
		E = 0.5*self.M*self.Pmass* (vx**2 + vy**2 + vz**2) / self.eCharge
		Ehit = E[hit]
		chargeHit = self.data['positions']['charge'][hit]
		
		for j in xrange(Ny):
			if j == Ny-1: fy = (ys >= yedges[j]) & (ys <= yedges[j+1])
			else: fy = (ys >= yedges[j]) & (ys < yedges[j+1])
			for i in xrange(Nx):
				if i == Nx-1: fx = (xs >= xedges[i]) & (xs <= xedges[i+1])
				else: fx = (xs >= xedges[i]) & (xs < xedges[i+1])
				fs = fx & fy
				energies[j,i] = Ehit[fs].mean()          # Note: mean of an empty array is nan
				charges[j,i] = chargeHit[fs].mean()

		fig = plt.figure(figsize = (10,6))		
		cs = plt.imshow(energies, extent = [X.min(), X.max(), Y.min(), Y.max()], cmap = 'plasma', origin = 'lower', aspect = 'auto', interpolation = 'nearest')
		plt.xlabel('x [m]')
		plt.ylabel('y [m]')
		C = plt.colorbar(cs, pad = 0.01, format = '%.3g')
		C.set_label('Particle Energy [eV]', rotation = 270, size = FONTSIZE, va = 'bottom')
		plt.title('Target Erosion: Energy', fontsize = FONTSIZE)

		fig = plt.figure(figsize = (10,6))		
		cs = plt.imshow(charges, extent = [X.min(), X.max(), Y.min(), Y.max()], cmap = 'plasma', origin = 'lower', aspect = 'auto', interpolation = 'nearest')
		plt.xlabel('x [m]')
		plt.ylabel('y [m]')
		C = plt.colorbar(cs, pad = 0.01, format = '%.3g')
		C.set_label('Particle Charge [e]', rotation = 270, size = FONTSIZE, va = 'bottom')
		plt.title('Target Erosion: Charge', fontsize = FONTSIZE)


	def plotOrbits(self, N = None):
		"""
		Plot the 3D orbits of the first N particles (default is all)
		"""
		if not self.data.has_key('history'): 
			print 'history.nc not available'
			return

		x = self.data['history']['x']
		y = self.data['history']['y']
		z = self.data['history']['z']
		r = np.sqrt(x**2 + y**2)
		
		fig = plt.figure(figsize = (9,6))
		ax = fig.add_subplot(111, projection='3d')
		if N is None: N = self.data['history']['nP']
		for i in xrange(N):
			ax.plot(x[i,:],y[i,:],z[i,:])
		plt.title(self.impurity + ' Particle Orbits', fontsize = FONTSIZE)
		ax.view_init(azim = 225)
		ax.ticklabel_format(style = 'sci',useOffset = True, useMathText = True, scilimits = (2,2))
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')
		
		plt.figure(figsize = (9,6))
		for i in xrange(N):
			plt.plot(r[i,:],z[i,:])
		plt.title(self.impurity + ' Particle Orbits', fontsize = FONTSIZE)
		plt.xlabel('r [m]')
		plt.ylabel('z [m]')


	def plotDensity(self):
		"""
		Plot the particle density in the X-section
		"""
		if not self.data.has_key('spec'): 
			print 'spec.nc not available'
			return

		dens =  self.data['spec']['n'][-1,:,:,:]
		gridR =  self.data['spec']['gridR']
		gridY =  self.data['spec']['gridY']
		gridZ =  self.data['spec']['gridZ']
		sumdens = np.sum(dens,1)/len(self.data['spec']['gridY'])
		digit = int(np.floor(np.log10(sumdens.max())))
		if abs(digit) > 1:
			factor = 10**digit
			sumdens /= factor
		sumdens[sumdens <= 0] = np.nan		# makes areas with zero density white in the plot
		fig = plt.figure(figsize = (10,6))		
		cs = plt.imshow(sumdens, extent = [gridR.min(), gridR.max(), gridZ.min(), gridZ.max()], cmap = 'plasma', origin = 'lower', aspect = 'auto', interpolation = 'nearest')
		plt.xlabel('x [m]')
		plt.ylabel('z [m]')
		C = plt.colorbar(cs, pad = 0.01, format = '%.3g')
		if abs(digit) > 1:
			C.set_label('Density [count x 10$^{' + str(digit) + '}$]', rotation = 270, size = FONTSIZE, va = 'bottom')
		else:
			C.set_label('Density [count]', rotation = 270, size = FONTSIZE, va = 'bottom')
		plt.title('X-section Particle Density', fontsize = FONTSIZE)
		plt.plot(gridR, gridR* np.sin(np.pi/6.0), 'k-', lw = 2)
		
	
	def plotEnergyDist(self, all = False):
		"""
		Plot the Energy-Angle distribution of the particles
		"""
		if not self.data.has_key('surface'): 
			print 'surface.nc not available'
			return

		if all:
			for key in self.data['surface'].keys():
				if key in ['nSurfaces','nAngles','nEnergies','surfEDist','surfaceNumber']: continue
				plt.figure(figsize = (9,6))
				plt.plot(self.data['surface']['surfaceNumber'],self.data['surface'][key], 'k-', lw = 2)
				plt.title(key, fontsize = FONTSIZE)
				plt.xlabel('Surface #')
				plt.ylabel('Value [a.u.]')
		
		edist = self.data['surface']['surfEDist'].flatten()
		edist = edist.reshape(self.data['surface']['nSurfaces'],self.data['surface']['nEnergies'],self.data['surface']['nAngles'])
		eDtotal = np.sum(edist,0).T
		#eDtotal *= 100.0/eDtotal.max()
		eDtotal[eDtotal <= 0] = np.nan		# makes areas with zero density white in the plot
		fig = plt.figure(figsize = (10,6))		
		cs = plt.imshow(eDtotal, extent = [0, 1000, 0, 90], cmap = 'plasma', origin = 'lower', aspect = 'auto', interpolation = 'nearest')
		plt.xlabel('Energy [eV]')
		plt.ylabel('Angle [deg]')
		C = plt.colorbar(cs, pad = 0.01, format = '%.3g')
		C.set_label('Particles [a.u.]', rotation = 270, size = FONTSIZE, va = 'bottom')
		plt.title('Energy Distribution', fontsize = FONTSIZE)
		

	def show(self, KEY = None):
		"""
		List the dimensions and variables in the file KEY. Default is: list all available files
		"""
		for file in self.files:
			if KEY is None: print '\n--------', file,'--------'
			else: file = KEY
			
			arrays,dims,other = [],[],[]
			for key in self.data[file].keys():
				if isinstance(self.data[file][key],np.ndarray): arrays.append(key)
				elif isinstance(self.data[file][key],int): dims.append(key)
				else: other.append(key)
			
			print 'Dimensions:'
			for key in dims: print key, '=', self.data[file][key]
		
			print '\n','Arrays:'
			for key in arrays: print key, self.data[file][key].shape

			if len(other) > 0:
				print '\n','Other Variables:'
				for key in other: print key, '=', self.data[file][key]
			
			if KEY is not None: break


	def modifyInputTimeSteps(self, nT = 100, path = None):
		"""
		Modify the number of time steps nT in the main input file
		"""
		if path is None: path  = self.path
		self.config['gitrInput']['timeStep']['nT'] = nT
		file = path + 'input/gitrInput.cfg'
		with io.open(file,'w') as f:
			libconf.dump(self.config['gitrInput'],f)


	def plot2dGeom(self, fig = None):
		"""
		Plot the 2D geometry
		"""
		if not self.config.has_key('gitrGeometry'): 
			print 'gitrGeometry.cfg not available'
			return

		x1 = np.array(self.config['gitrGeometry'].geom.x1)	# both notations, with . or with [''] are possible
		x2 = np.array(self.config['gitrGeometry'].geom.x2)
		z1 = np.array(self.config['gitrGeometry'].geom.z1)
		z2 = np.array(self.config['gitrGeometry'].geom.z2)
		Z = np.array(self.config['gitrGeometry'].geom.Z)
		length = np.array(self.config['gitrGeometry'].geom.length)
		
		if fig is None: 
			plt.figure(figsize = (9,6))
			plt.title('Geometry', fontsize = FONTSIZE)
			plt.xlabel('x [m]')
			plt.ylabel('z [m]')
		plt.plot(np.append(x1,x1[0]),np.append(z1,z1[0]),'k-',lw = 2)
		plt.xlim(x1.min()-0.05*abs(x1.min()),x1.max()+0.05*abs(x1.max()))
		plt.ylim(z1.min()-0.05*abs(z1.min()),z1.max()+0.05*abs(z1.max()))


	def plot3dGeom(self):
		"""
		Plot the 3D geometry
		"""
		if not self.config.has_key('gitrGeometry'): 
			print 'gitrGeometry.cfg not available'
			return

		x1 = np.array(self.config['gitrGeometry'].geom.x1)
		x2 = np.array(self.config['gitrGeometry'].geom.x2)
		x3 = np.array(self.config['gitrGeometry'].geom.x3)
		y1 = np.array(self.config['gitrGeometry'].geom.y1)
		y2 = np.array(self.config['gitrGeometry'].geom.y2)
		y3 = np.array(self.config['gitrGeometry'].geom.y3)
		z1 = np.array(self.config['gitrGeometry'].geom.z1)
		z2 = np.array(self.config['gitrGeometry'].geom.z2)
		z3 = np.array(self.config['gitrGeometry'].geom.z3)
		area = np.array(self.config['gitrGeometry'].geom.area)
		surf = np.array(self.config['gitrGeometry'].geom.surface)
		Z = np.array(self.config['gitrGeometry'].geom.Z)

		xs=[]
		ys=[]
		zs=[]
		for i in range(0,x1.size - 1):
			xs.append(x1[i])
			xs.append(x2[i])
			xs.append(x3[i])
			ys.append(y1[i])
			ys.append(y2[i])
			ys.append(y3[i])
			zs.append(z1[i])
			zs.append(z2[i])
			zs.append(z3[i])
		
		verts = [zip(xs,ys,zs)]

		fig = plt.figure(figsize = (9,6))
		ax = fig.add_subplot(111, projection='3d')
		ax.add_collection3d(Poly3DCollection(verts))
		plt.title('3D Geometry', fontsize = FONTSIZE)
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')
		ax.set_xlim3d(min(xs)-0.25*abs(min(xs)),max(xs)+0.25*abs(max(xs)))
		ax.set_ylim3d(min(ys)-0.25*abs(min(ys)),max(ys)+0.25*abs(max(ys)))
		ax.set_zlim3d(min(zs)-0.25*abs(min(zs)),max(zs)+0.25*abs(max(zs)))
		
		materialSurfaceInidces = np.nonzero(Z)
		surfIndArray = np.asarray(materialSurfaceInidces)
		print 'Number of W surfaces', surfIndArray.size
		

	def plotVz(self, N = None):
		"""
		Plot the vertical velocity of N orbits (default is all)
		"""
		if not self.data.has_key('history'): 
			print 'history.nc not available'
			return

		z = self.data['history']['z']
		vz = self.data['history']['vz']

		fig = plt.figure(figsize = (9,6))
		if N is None: N = self.data['history']['nP']
		for i in xrange(N):
			plt.plot(z[i,:],vz[i,:])

		plt.xlabel("z [m]")
		plt.ylabel("v$_z$ [m/s]")
		

	def plotPitch(self):
		"""
		Plot the pitch angle distribution and velocity distributions among the particles
		"""
		if not self.data.has_key('positions'): 
			print 'positions.nc not available'
			return

		vx = self.data['positions']['vx']
		vy = self.data['positions']['vy']
		vz = self.data['positions']['vz']
		nP = self.data['positions']['nP']
		weights = np.ones(nP)/np.float64(nP)*100
		
		vperp = np.sqrt(vx**2 + vy**2)
		pitchAngle = np.arctan(vperp/vz)
			
		fig = plt.figure(figsize = (9,6))		
		plt.hist(pitchAngle, bins = 30, weights = weights, rwidth = 0.8)
		plt.xlabel("Pitch Angle [rad]")
		plt.ylabel("# Particles [%]")
		plt.xlim(-1.6,1.6)
		
		fig = plt.figure(figsize = (9,12))
		plt.subplot(3,1,1)
		plt.hist(vx, bins = 30, weights = weights, rwidth = 0.8)
		plt.xlabel("v$_x$ [m/s]")
		plt.ylabel("# Particles [%]")
		
		plt.subplot(3,1,2)
		plt.hist(vy, bins = 30, weights = weights, rwidth = 0.8)
		plt.xlabel("v$_y$ [m/s]")
		plt.ylabel("# Particles [%]")
		
		plt.subplot(3,1,3)
		plt.hist(vz, bins = 30, weights = weights, rwidth = 0.8)
		plt.xlabel("v$_z$ [m/s]")
		plt.ylabel("# Particles [%]")
		
		
	def plotErosion(self, logplot = True):
		"""
		Plot the Deposition and Erosion along the surface
		"""
		if not self.data.has_key('surface'): 
			print 'surface.nc not available'
			return

		length = np.array(self.config['gitrGeometry'].geom.length)
		N = self.data['surface']['nSurfaces']
		surfaces = self.data['surface']['surfaceNumber']*length[0]/np.float64(N)
		grossEro = self.data['surface']['grossErosion']
		grossDep = self.data['surface']['grossDeposition']
		netErosion = grossEro - grossDep
		positiv = netErosion >= 0
		
		top = max([grossEro.max(),grossDep.max(),np.abs(netErosion).max()])
		digit = int(np.ceil(np.log10(top)))
			
		fig = plt.figure(figsize = (9,6))
		if logplot:
			s1 = plt.semilogy(surfaces, grossEro, 'bo', label = 'Gross Erosion')
			s2 = plt.semilogy(surfaces, grossDep, 'ro', label = 'Gross Deposition')
			s3 = plt.semilogy(surfaces[positiv], netErosion[positiv], 'go', label = 'net Erosion')
			s4 = plt.semilogy(surfaces[-positiv], -netErosion[-positiv], 'g^', label = 'net Deposition')
			plt.ylim(1e-2,10**digit)
		else:
			s1 = plt.plot(surfaces, grossEro, 'bo', label = 'Gross Erosion')
			s2 = plt.plot(surfaces, grossDep, 'ro', label = 'Gross Deposition')
			s3 = plt.plot(surfaces, netErosion, 'go', label = 'net Erosion')
			
		plt.legend(prop={'size': 12})
		plt.xlabel("s [m]")
		plt.ylabel('Rate [a.u.]')
		
		
	def plotChargeHistogram(self):
		"""
		Plot the charge distribution among the particles
		"""
		if not self.data.has_key('positions'): 
			print 'positions.nc not available'
			return
			
		fig = plt.figure(figsize = (9,6))
		bins = np.arange(self.data['positions']['charge'].max()+2)-0.5
		nP = self.data['positions']['nP']
		weights = np.ones(nP)/np.float64(nP)*100
		plt.hist(self.data['positions']['charge'],bins = bins, weights = weights, rwidth = 0.8)
		plt.xlabel('Charge [e]')
		plt.ylabel('# Particles [%]')

		

# ----------------------------------------------------------------------------------------
# --- End of Class -----------------------------------------------------------------------

# --- openNCfile(Filename) ---------------------------------------------------------------

def openNCfile(Filename):
	"""
	Open netcdf file 'Filename' and parse its contents into dictionary
	return dictionary
	"""
	from netCDF4 import Dataset
	from warnings import filterwarnings
	filterwarnings('ignore')
	
	D = {}
	data = Dataset(Filename, 'r')
	
	# Dimensions
	varnames = data.dimensions.keys()
	for n in varnames:
		D[n] = len(data.dimensions[n])
	
	# Variables
	varnames = data.variables.keys()
	for n in varnames:
		if(data.variables[n].size > 1):										# arrays (all types)
			D[n] = np.array(data.variables[n][:])
			if(data.variables[n].dtype == 'S1'):							# character array
				if(data.variables[n].size == data.variables[n].shape[0]):	# 1-D						
					D[n] = ''.join(data.variables[n][:]).strip()
				else:														# 2-D
					D[n] = []
					for i in xrange(data.variables[n].shape[0]):
						D[n].append(''.join(data.variables[n][i,:]).strip())
		else:																# single variable
			if(data.variables[n].dtype in ['float','float32','float64']):	# float
				D[n] = np.float64(data.variables[n][:])
			elif(data.variables[n].dtype in ['int','int32','int64','long']):# int
				D[n] = int(data.variables[n][:])
			elif(data.variables[n].dtype == 'S1'):							# char
				try: D[n] = ''.join(data.variables[n][:])	# at fixed bndy: mgrid_mode exists but is masked -> error
				except: D[n] = data.variables[n][:]
			else:															# unknown
				print 'unknown datatype in', Filename, ', variable:', n, ' type:',data.variables[n].dtype, 'size:', data.variables[n].size
				D[n] = data.variables[n][:]
	return D


# -------------------------------------------------------------------------------------------------------------
# --- main ----------------------------------------------------------------------------------------------------
def main():
	g = gitr(show = True)
	plt.show()

if __name__ == '__main__':
	import argparse
	import textwrap
	parser = argparse.ArgumentParser(description = 'Load, Post-process and Plot GITR output', 
				formatter_class = argparse.RawDescriptionHelpFormatter,
				epilog = 'Please report bugs to: wingen@fusion.gat.com')

	parser.add_argument('-v', '--version', help = 'Show current version number and release data', action = 'store_true', default = False)
	args = parser.parse_args()

	if args.version: print '\ngitr_class, Version: ' + str(_VERSION) + ', Release: ' + _LAST_UPDATE + '\n'
	else: main()



# ----------------------------------------------------------------------------------------
# --- Old legacy code --------------------------------------------------------------------


"""
#from distutils.dir_util import copy_tree
#import sys
#sys.path.append('/home/tqd/code/netcdf4-python')

#import Tkinter

#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.use('agg')
#import cv2
#import io,libconf

#import pylab as pl
#import scipy as sp
#import math
#import shutil


    

def iter2dProcessing():
    x1,x2,z1,z2,length,Z = plot2dGeom('input/iterRefinedTest.cfg')
    plt.close()
    grossDep,grossEro,sumWeightStrike,E,A,EAdist,surfaceNumbers,surfEdist = nc_readSurface()
    plt.close()
    netErosion = grossDep - grossEro
    colormap = plt.cm.bwr
    normalize = matplotlib.colors.Normalize(vmin=-1.,vmax=1.) 
    plt.subplot(2,1,1)
    for i in range(0,len(surfaceNumbers)):
        plt.plot([x1[surfaceNumbers[i]], x2[surfaceNumbers[i]]], [z1[surfaceNumbers[i]], z2[surfaceNumbers[i]]],color=colormap(0.))
        plt.hold(True)
    plt.subplot(2,1,2)
    plt.scatter(surfaceNumbers,netErosion)
    plt.savefig('surfaces.png') 
    plt.close()
    x,y,r,z,charge = nc_plotPositions('output/positions.nc')
    plt.hist(charge)
    plt.savefig('charge.png')
    plt.close
    rSep = 5.5543
    zSep = -4.3970
    gridRsep = np.zeros(160)
    gridZsep = np.zeros(160)
    rMinusRsep = np.zeros(160)
    for i in range(160):
        print(i)
        thisSurf = surfaceNumbers[i+4]
        lastSurf = surfaceNumbers[i+3]
        print(thisSurf)
        print(z2[thisSurf])
        if(i == 0):
            rMinusRsep[i] = z2[thisSurf] - zSep
        else:
            rMinusRsep[i] = rMinusRsep[i-1] + length[lastSurf]
    plt.close('all')
    s1 = plt.scatter(rMinusRsep,grossEro[range(4,164)],c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep,grossDep[range(4,164)],c='red')
    s3 = plt.scatter(rMinusRsep,netErosion[range(4,164)],c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosion.png') 
    plt.close('all')
    s1 = plt.scatter(rMinusRsep,np.log10(grossEro[range(4,164)]),c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep,np.log10(grossDep[range(4,164)]),c='red')
    s3 = plt.scatter(rMinusRsep,np.sign(netErosion[range(4,164)])*np.log10(np.absolute(netErosion[range(4,164)])),c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosionlog.png') 
    plt.close()
    
def printHeDist(path = '',z = [-4.1,-4.0]):
    ncFile = netCDF4.Dataset("input/iterHeDist.nc","r")
    nPoints = len(ncFile.dimensions['nPoints'])
    nE = len(ncFile.dimensions['nE'])
    nA = len(ncFile.dimensions['nA'])

    zPoints = np.array(ncFile.variables['z'])
    flux = np.array(ncFile.variables['flux'])
    EAdist = np.array(ncFile.variables['EAdist'])
    gridE = np.array(ncFile.variables['gridE'])
    gridA = np.array(ncFile.variables['gridA'])
    heFlux = np.zeros(len(z))
    os.mkdir('GITRoutput')
    for i in range(len(z)):
        idx = (np.abs(zPoints - z[i])).argmin()
        print('index',idx)
        heFlux[i] = flux[idx] 
        thisDist = EAdist[:,:,idx] #first dimension is angle, second is energy
        thisDist = np.transpose(thisDist)
        print('thisDistHe size', thisDist.shape)
	#print(thisDist)
        Aweight = np.sum(thisDist,axis=0)
        gitrDir = 'GITRoutput/'+'gitr'+str(i)
        os.mkdir(gitrDir)
        np.savetxt(gitrDir+'/gitrFluxE.dat', gridE)
        np.savetxt(gitrDir+'/gitrFluxAweight.dat', Aweight)
        np.savetxt(gitrDir+'/gitrFluxA.dat',gridA[:-1])
        np.savetxt(gitrDir+'/gitrFluxEAdist.dat', thisDist)
        
        if(path != ''): 
            np.savetxt(path+'/'+gitrDir+'/gitrFluxE.dat', gridE)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxAweight.dat', Aweight)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxA.dat', gridA)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxEAdist.dat', thisDist)

        for j in range(0,nA):
            #print('Evec size ', gridE.size)
            #print('EAdist[:,i] size ', EAdist[:,i].size)
            edOut = np.column_stack((gridE,thisDist[:,j]))
            np.savetxt(gitrDir+'/dist'+str(j)+'.dat', edOut)
    return gitrDir, heFlux
    
def iter3dProcessing(path = '',loc = [-4.5020,   -4.4628,   -4.4237,   -4.3846,   -4.3455,   -4.3064,   -4.2672,   -4.2281, -4.1890,   -4.1499,   -4.1108,   -4.0717,   -4.0327,   -3.9938,   -3.9547,   -3.9157, -3.8766,   -3.8380,   -3.7395,   -3.6386,   -3.5423,   -3.4537,   -3.3799,   -3.3195,-3.2744],locWidth = 0.02):
    #loc = [-4.5020,   -4.4628,   -4.4237,   -4.3846,   -4.3455,   -4.3064,   -4.2672,   -4.2281, -4.1890,   -4.1499,   -4.1108,   -4.0717,   -4.0327,   -3.9938   -3.9547,   -3.9157, -3.8766,   -3.8380,   -3.7395,   -3.6386,   -3.5423,   -3.4537,   -3.3799,   -3.3195,-3.2744]
    #loc = [-4.50,-4.437,-4.3570,-4.2770,-4.1970,-4.1170,-4.0371,-3.9574,-3.8775,-3.5341,-3.3679,-3.273]
    #x1,x2,z1,z2,length,Z = plot2dGeom('input/iterRefinedTest.cfg')
    x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray,surf = read3dGeom('input/iterRefinedTest.cfg')
    plt.close()
    grossDep,grossEro,sumWeightStrike,E,A,EAdist,surfaceNumbers,surfEdist = nc_readSurface()
    plt.close()
    print('e', E)
    #loc1 = -4.17
    #loc2 = -4.081
    #loc3 = -4.25
    #loc4 = -3.6
    #locWidth = 0.02
    surfInd = surf >0
    z1= np.extract(surfInd,z1)
    gitrDirHe,heFlux=printHeDist(z=loc)
    os.mkdir('GITRoutput_W')
    nLocations = len(loc)
    for i in range(len(loc)):
        condition = [(z1 < loc[i]+locWidth) & (z1 > loc[i]-locWidth)]
	theseInd = np.where(condition)[1]
	print('theseInd',theseInd)
	print('condition',condition)
	print('surfEdist Shape',surfEdist.shape)
	print('theseInd',theseInd)
	print('theseInd',theseInd.shape)
	thisDist = surfEdist[theseInd,:,:]
	print('size thisDist',thisDist.shape)
	thisDist = np.sum(thisDist,axis=0)
	print('size thisDist',thisDist.shape) 
	#thisDist = np.sum(thisDist,axis=1)
	print('size thisDist',thisDist.shape) 
	print('thisDist',thisDist) 
        dep = np.extract(condition,grossDep)
        ero = np.extract(condition,grossEro)
        strike = np.extract(condition,sumWeightStrike)
        areas = np.extract(condition,area)
        with io.open('input/gitrInput.cfg') as f:
            config = libconf.load(f)
        #backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
        #backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
        #time = float(config.postProcessing.time);
        nParticles = float(config.impurityParticleSource.nP);
        #backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
        erodedFlux = float(config.postProcessing.totalWFlux);
        erodedFluxPerParticle = erodedFlux/nParticles;
        netErosion = np.sum(ero - dep);
        netStrike = np.sum(strike)
        totalArea = np.sum(areas)
        impurityFlux = netErosion*erodedFluxPerParticle;
        if(heFlux[i]==0.0):
            heFlux[i] = 1.0e19
        Wfrac = impurityFlux/heFlux[i];
        Aweight = np.sum(thisDist,axis=1)
        print('W impurity flux ', impurityFlux)
        print('W impurity fraction ', Wfrac)
        #for i in surfIndArray:
	gitrDir = 'GITRoutput_W/'+'gitr'+str(i)
	os.mkdir(gitrDir)
        np.savetxt(gitrDir+'/gitrFluxE.dat', E)
        np.savetxt(gitrDir+'/gitrFluxAweight.dat', Aweight)
        np.savetxt(gitrDir+'/gitrFluxA.dat', A[:-1])
        np.savetxt(gitrDir+'/gitrFluxEAdist.dat', thisDist)
        
        if(path != ''): 
            np.savetxt(path+'/'+gitrDir+'/gitrFluxE.dat', E)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxAweight.dat', Aweight)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxA.dat', A)
            np.savetxt(path+'/'+gitrDir+'/gitrFluxEAdist.dat', thisDist)
        
        Dfrac = float(config.postProcessing.Dfrac);
        Hefrac = float(config.postProcessing.Hefrac);
        Tfrac = float(config.postProcessing.Tfrac);
        file = open('gitrOut_'+str(i)+'.txt','w') 
        file.write('plasmaSpecies=He W D T\n') 
	file.write('inputEnergy=-1.0 -1.0 0.0 0.0\n')
	file.write('inputAngle=-1.0 -1.0 0.0 0.0\n')
        file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
        file.write('flux='+str(heFlux[i]/1e18)+'\n') 
        file.write('gitrOutputDir_He='+os.getcwd()+'/'+'GITRoutput/'+'gitr'+str(i)+'\n') 
        file.write('gitrOutputDir_W='+os.getcwd()+'/'+gitrDir+'\n') 
        file.close() 

        if(path != ''): 
            shutil.copyfile('gitrOut.txt',path+'/'+gitrDir+'/gitrOut.txt')
            #file = open(path+'/'+'gitrOut.txt','w') 
            #
            #file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
            #file.write('flux='+str(backgroundFlux/1e18)+'\n') 
            #file.write('gitrOutputDir='+os.getcwd()+'\n') 
            #file.close() 

        E0=0.0;
        E1=1000.0;
        nE=200;
        dE = E1/nE;
        Evec = np.linspace(0.5*dE,E1-0.5*dE,nE);
        A0=0.0;
        A1=90.0;
        nA=30;
        for j in range(0,nA):
            #print('Evec size ', Evec.size)
            #print('EAdist[:,i] size ', EAdist[:,i].size)
            edOut = np.column_stack((Evec,EAdist[:,j]))
            np.savetxt(gitrDir+'/dist'+str(j)+'.dat', edOut)
    return nLocations
    
def piscesProcessing(r=0.01,path=''):
    x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray = read3dGeom('input/gitrGeometryPisces1inch.cfg')
    r1 = np.sqrt(np.multiply(x1[surfIndArray],x1[surfIndArray]) + np.multiply(y1[surfIndArray],y1[surfIndArray]))
    r2 = np.sqrt(np.multiply(x2[surfIndArray],x2[surfIndArray]) + np.multiply(y2[surfIndArray],y2[surfIndArray]))
    r3 = np.sqrt(np.multiply(x3[surfIndArray],x3[surfIndArray]) + np.multiply(y3[surfIndArray],y3[surfIndArray]))
    grossDep,grossEro,sumWeightStrike,E,A,EAdist = nc_readSurface()
    condition = r1 < r
    dep = np.extract(condition,grossDep)
    ero = np.extract(condition,grossEro)
    strike = np.extract(condition,sumWeightStrike)
    areas = np.extract(condition,area)
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)

    backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
    backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
    time = float(config.postProcessing.time);
    nParticles = float(config.impurityParticleSource.nP);
    backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    erodedMass = time*backgroundIonsPerSec*184*1.66e-27*backgroundSputtYield*1000;
    erodedMassPerParticle = erodedMass/nParticles;
    netErosion = np.sum(ero - dep);
    netStrike = np.sum(strike)
    totalArea = np.sum(areas)
    impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec*backgroundSputtYield/nParticles;
    impurityFlux = netStrike/totalArea*impurityParticlePerSecondPerComputationalPartice;
    Wfrac = impurityFlux/backgroundFlux;
    Aweight = np.sum(EAdist,axis=0)
    print('W impurity flux ', impurityFlux)
    print('W impurity fraction ', Wfrac)
    #for i in surfIndArray:
    np.savetxt('gitrFluxE.dat', E)
    np.savetxt('gitrFluxAweight.dat', Aweight)
    np.savetxt('gitrFluxA.dat', A[:-1])
    np.savetxt('gitrFluxEAdist.dat', EAdist)
    
    if(path != ''): 
        np.savetxt(path+'/'+'gitrFluxE.dat', E)
        np.savetxt(path+'/'+'gitrFluxAweight.dat', Aweight)
        np.savetxt(path+'/'+'gitrFluxA.dat', A)
        np.savetxt(path+'/'+'gitrFluxEAdist.dat', EAdist)
    
    Dfrac = float(config.postProcessing.Dfrac);
    Hefrac = float(config.postProcessing.Hefrac);
    Tfrac = float(config.postProcessing.Tfrac);
    file = open('gitrOut.txt','w') 
    file.write('plasmaSpecies=He W D T\n') 
    file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
    file.write('flux='+str(backgroundFlux/1e18)+'\n') 
    file.write('gitrOutputDir='+os.getcwd()+'\n') 
    file.close() 

    if(path != ''): 
        shutil.copyfile('gitrOut.txt',path+'/'+'gitrOut.txt')
        #file = open(path+'/'+'gitrOut.txt','w') 
        #
        #file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
        #file.write('flux='+str(backgroundFlux/1e18)+'\n') 
        #file.write('gitrOutputDir='+os.getcwd()+'\n') 
        #file.close() 

    E0=0.0;
    E=1000.0;
    nE=100;
    dE = E/nE;
    Evec = np.linspace(0.5*dE,E-0.5*dE,nE);
    A0=0.0;
    A=90.0;
    nA=90;
    for i in range(0,nA):
        print('Evec size ', Evec.size)
        print('EAdist[:,i] size ', EAdist[:,i].size)
        edOut = np.column_stack((Evec,EAdist[:,i]))
        np.savetxt('dist'+str(i)+'.dat', edOut)

def d3dProcessing(r=0.01,path=''):
    #plt.figure,figsize=(6, 10)) #, dpi=1000)
    #plot2dGeom(filename='../input/gitrD3DGeometry2DWrings.cfg')
    grossDep,grossEro,sumWeightStrike,E,A,EAdist = nc_readSurface()
    print('W gross deposition ', grossDep)
    print('W gross erosion ', grossEro)
    #condition = r1 < r
    #dep = np.extract(condition,grossDep)
    #ero = np.extract(condition,grossEro)
    #strike = np.extract(condition,sumWeightStrike)
    #areas = np.extract(condition,area)
    #with io.open('input/gitrInput.cfg') as f:
    #    config = libconf.load(f)

    #backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
    #backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
    #time = float(config.postProcessing.time);
    #nParticles = float(config.impurityParticleSource.nP);
    #backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    #erodedMass = time*backgroundIonsPerSec*184*1.66e-27*backgroundSputtYield*1000;
    #erodedMassPerParticle = erodedMass/nParticles;
    #netErosion = np.sum(ero - dep);
    #netStrike = np.sum(strike)
    #totalArea = np.sum(areas)
    #impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec*backgroundSputtYield/nParticles;
    #impurityFlux = netStrike/totalArea*impurityParticlePerSecondPerComputationalPartice;
    #Wfrac = impurityFlux/backgroundFlux;
    #Aweight = np.sum(EAdist,axis=0)
    #print('W impurity flux ', impurityFlux)
    #print('W impurity fraction ', Wfrac)





if __name__ == "__main__":
    #asdfanc_show("surface.nc")
    #depositedEdist()
    #if(os.path.exists('output/history.nc')):
    # 	nc_plotHist('output/history.nc')
    #if(os.path.exists('output/spec.nc')):
    #	nc_plotSpec('output/spec.nc')
    #iter2dProcessing()
    iter3dProcessing()
    #printHeDist()
    #nc_plotSpec3D()
    #nc_plotPositions()
    #nc_plotVz()
    #plotPitch()
    #piscesProcessing()
    #modifyInputParam()
    #nc_readSurface()
"""