"""
Tim Waters
waterst3@unlv.nevada.edu

--------- DESCRIPTION ---------
This class is a reader for 1D Athena4.2 tab output.
It can be used to open the default output files, as
well as the user-defined '.out' files.  
"""

import numpy as np
import sys
		
class AthenaSim1D:
	def __init__(self,fileinfo):
		#unpack fileinfo dictionary
		self.simpath = fileinfo['simpath']
		self.foldername = fileinfo['simfolder']
		self.simtag = fileinfo['simtag']
		self.g = fileinfo['gamma']
		
		# list to hold time extension in between updates
		self.timestamp = ['','','','']
		
	def restart_check(self,lev):
		restart = self.timestamp[lev]
		smrtag = ''
		if lev != 0:
			smrtag = '-lev' + `lev`
		name = '/' + self.simtag + smrtag
		if self.i != 0 and self.i%(Nrst+1) == 0:
			self.i=Nrst
			if self.i < 10:
				restart = name + '.000'
			elif self.i < 100:
				restart = name + '.00'
			elif self.i < 1000:
				restart = name + '.0'
			elif self.i < 10000:
				restart = name + '.'
		return restart
		
	def update_check(self,lev):
		update = self.timestamp[lev]
		smrtag = ''
		if lev != 0:
			smrtag = '-lev' + `lev`
		name = '/' + self.simtag + smrtag
		if self.i < 10:
			update = name + '.000'
		if self.i>=10: 
			update = name + '.00'
		if self.i>=100:
			update = name + '.0'
		if self.i>=1000:
			update = name + '.'
		return update
		
	def get_filepath(self,i,ext,lev=0,Nrst=10000):
		self.i = i
		smrtag = ''
		if lev != 0:
			smrtag = '/lev' + `lev`
		self.timestamp[lev] = self.update_check(lev)
		if Nrst > 0 and Nrst < 10000:
			self.timestamp[lev] = self.restart_check(lev)
		return self.simpath + self.foldername + smrtag + self.timestamp[lev] + `self.i` + ext
		
	def load_data(self,i,lev=0,out_file='',start=0,end=0,Nrst=10000):
		 # to load a custom athena .out file, pass an integer to out_file
		 # e.g. to load _.out2.tab pass out_file = 2
		if out_file != '':
			ext = out_file + '.tab'
		else:
			ext = '.tab'
		filepath = self.get_filepath(i,ext,lev,Nrst)
		self.data = np.loadtxt(filepath)
		
	def hydro_profile(self,var,start=0,end=0):
		if var == 'index':
			y = self.data[:,0]	
		elif var == 'x':
			y = self.data[:,1]
		elif var == 'density':
			y = self.data[:,2]
		elif var == 'velocity':
			y = self.data[:,3]
		elif var == 'pressure':
			y = self.data[:,6]
		
		# derived quantities
		elif var == 'momentum':
			rho = self.data[:,2]
			v = self.data[:,3]
			y = rho*v
		elif var == 'entropy':
			P = self.data[:,6]
			rho = self.data[:,2]
			y = np.log(P/rho**self.g)
		elif var == 'sound speed':
			P = self.data[:,6]
			rho = self.data[:,2]
			y = np.sqrt(self.g*P/rho)
		elif var == 'Mach number':
			v = self.data[:,3]
			P = self.data[:,6]
			rho = self.data[:,2]
			cs = np.sqrt(self.g*P/rho)
			y = v/cs
		elif var == 'internal energy' and self.g > 1.0:
			P = self.data[:,6]
			y = P/(self.g-1.)	
		elif var == 'total energy' and self.g > 1.0:
			v = self.data[:,3]
			P = self.data[:,6]
			rho = self.data[:,2]
			e = P/(self.g-1.)
			y = e + 0.5*rho*v**2	
		else:
			print "\nERROR in AthenaSim1D!"
			print "Didn't recognize variable = {}".format(var)
			sys.exit()
			
		if end == 0:
			end = len(y)
		return y[start:end]	