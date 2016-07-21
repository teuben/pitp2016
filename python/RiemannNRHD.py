"""
Timothy R. Waters <waterst3@unlv.nevada.edu>
University of Nevada Las Vegas, Summer 2013
(GNU General Public License)


--------- PROGRAM SUMMARY ---------
The 'NRHD' in RiemannNRHD.py stands for 'non-relativistic hydrodynamic', as the equations
solved, namely the 1D Euler equations, are in this physical regime. RiemannNRHD solves the 
1D Riemann problem based on the solution found in chapter 4 of
	E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics,
	DOI 10.1007/b7976-1_1 Springer-Verlag Berlin Heidelberg 2009
All references to "Toro" refer to this book.  This solution assumes an ideal gas equation
of state (EOS), in which Internal_Energy = Pressure/((Gamma - 1)*density).  The isothermal
case, in which Gamma = 1 and Internal_Energy = Infinity, is also solved following the 
methods of Toro, although the actual equations are not given in that book.  (For the most
part, any terms involving 2/(Gamma - 1) were replaced by ln(density).)

*** Class RiemannProblem(left_state,right_state) ***
This class computes the solution to the Riemann problem for a given set of left and right
ambient states.  A state is defined by its adiabatic index, density, velocity, and 
pressure OR sound speed.  For example, the benchmark Sod shock tube is solved as follows:

WL = State(1.4,1.0,0.0,1.0)   #gamma=1.4, density=1.0,   velocity=0., pressure=1.0
WR = State(1.4,0.125,0.0,0.1) #gamma=1.4, density=0.125, velocity=0., pressure=0.1
sod = RiemannProblem(WL,WR)
Solution = sod.solution()

To run other Riemann problems, simply change the input states WL and WR. 

*** Class Exact1DProfile(Solution,xaxis,time,center=0.,resize=True) ***
This class requires the output (Solution) of the solution method of class
RiemannProblem and a given range of x values (xaxis) at given fixed moment in time (time).
The location of the interface between the left and right states is taken to be 0 by default
but can be changed using the optional argument "center".  The optional argument "resize" 
can be set to False to prevent adding the precise locations of the wavefronts to xaxis.  
Appending the precise locations of discontinuities makes the profiles look mathematically 
exact, but this is undesired if, say, doing a convergence test and the exact solution is 
desired at only the points contained in the original xaxis input array.  

By calling Exact1DProfile at many different time steps, the time-dependent evolution of
any available hydrodynamical variable can be visualized and used to compare with numerical
simulations...an example animation is provided at the bottom of this script. 

--------- ILLUSTRATIVE EXAMPLES ---------
Other interesting sets of states are
(1) Vacuum left state with zero density and pressure combined with a 
non-vacuum nearly isothermal right state.
WL = (5./3.,0.,1.,0.) 
WL = (1.01,0.1,0.,0.1) 

(2) Pressure-positivity violating states that have no 'star-region' and generate a vacuum
    region in between two rarefaction waves.
WL = State(5./3.,0.1,0.5,a=0.5) #sound speed is specified instead of pressure
WR = State(5./3.,0.125,10.,a=0.1) 

(3) Colliding pressure equilibrium gases of different adiabatic indices
Case 1. Double Shocks
WL = State(5./3.,1.,2.,1.)
WR = State(1.,1.,0.,1.)
Case 2. Double Rarefactions
WL = State(5./3.,1.,0.,1.)
WR = State(1.,1.,2.,1.)
Case 3 Flatline: If no relative velocity, 
There will be no contact discontinuity in each case if the adiabatic index is the same.

"""
import sys
import numpy as np
import pylab as plt
from scipy.optimize import brentq

# Struct to hold a scalar state of density, velocity, & pressure
class State:
	def __init__(self,gamma,d,u,p=None,a=None):
		self.d = d
		self.u = u
		self.g = gamma
		if a == None:
			if p == None:
				print "\nINVALID STATE: Must specify pressure or sound speed!"
				sys.exit()
			else:
				self.p = p 
				if d == 0.: self.a = 0.
				else: self.a = np.sqrt(self.g*self.p/self.d)
		if p == None:
			if a == None:
				print "\nINVALID STATE: Must specify pressure or sound speed!"
				sys.exit()
			else:
				self.a = a
				self.p = (a**2.)*self.d/self.g
		
# All methods required to completely solve the Riemann problem for an ideal gas
class RiemannProblem:	
	def __init__(self,WL,WR):
		print  "\n*** Inputs to Riemann Problem: ***"
		print "Left State:  (gamma = {}, density = {}, velocity = {}, pressure = {}, sound speed = {})"\
				.format(WL.g,WL.d,WL.u,WL.p,WL.a)
		print "Right State: (gamma = {}, density = {}, velocity = {}, pressure = {}, sound speed = {})"\
				.format(WR.g,WR.d,WR.u,WR.p,WR.a)
		if WL.g == 1.0 and WR.g == 1.0 and WL.a != WR.a:
			print '\nERROR: Invalid isothermal input states.' \
					'\nLeft and right sound speeds not equal!' 
			sys.exit()
		if WL.d == 0.0 and WR.d == 0.0 or WL.p == 0.0 and WR.p == 0.0:
			print '\nERROR: One state must be non-vacuum!' 
			sys.exit()
		if WL.d == 0.0 and WL.p != 0.0 or WR.d == 0.0 and WR.p != 0.0:
			print '\nERROR: Density = 0 implies pressure = 0 for an ideal gas!' 
			sys.exit()
		if WL.p == 0.0 and WL.d != 0.0 or WR.p == 0.0 and WR.d != 0.0:
			print '\nERROR: Pressure = 0 implies density = 0 for an ideal gas!' 
			sys.exit()
		if WL.p == 0.0 and WL.d == 0.0 and WR.g == 1.0:
			print '\nERROR: No solution for vacuum opposite an isothermal state!' 
			sys.exit()
		if WR.p == 0.0 and WR.d == 0.0 and WL.g == 1.0:
			print '\nERROR: No solution for vacuum opposite an isothermal state!'
			sys.exit()
		self.WL,self.WR = WL,WR
		self.u_star = None
		self.p_star = None
		self.dL_star = None
		self.dR_star = None
		
	### Check for vacuum states or a generated vacuum region ###
	def vacuum(self):
		#if self.WL.g == 1. or self.WR.g == 1.: return False #no vacuum if isothermal!
		if self.WL.p == 0. or self.WR.p == 0.: return True
		if self.WL.d == 0. or self.WR.d == 0.: return True
		if self.WL.g == 1. or self.WR.g == 1.: return False #no vacuum if isothermal!
		if self.WL.g == 1.0:
			facL = np.log(self.WL.d)
		else:
			facL = 2./(self.WL.g - 1.)
		if self.WR.g == 1.0:
			facR = np.log(self.WR.d)
		else:
			facR = 2./(self.WR.g - 1.)
		udelta = self.WR.u - self.WL.u
		pressureFunc = facL*self.WL.a + facR*self.WR.a - udelta
		if pressureFunc < 0.: 
			print "\n*** Pressure positivity condition violated! ***"
			print "Vacuum region generated in between rarefaction waves."
			print "Pressure function = {} (see Toro eqn. [4.40])".format(pressureFunc)
			return True
		else: return False
	
	### Build Pressure Function: Next 5 functions define (4.5) in Toro ###
	def f_shock(self,W,p):
		gm1 = W.g - 1.
		gp1_inv = 1./(W.g + 1.)
		aK = 2.*gp1_inv/W.d
		bK = W.p*gm1*gp1_inv
		return (p-W.p)*np.sqrt(aK/(p+bK))
		
	def f_rarefaction(self,W,p):
		if W.g == 1.0:
			return W.a*(np.log(p) - np.log(W.p))
		else:
			gm1 = W.g - 1.
			brak = (p/W.p)**(0.5*gm1/W.g) - 1.
			return 2.*W.a*brak/gm1
		
	def fL(self,p):
		if p>self.WL.p and self.WL.p != 0.: return self.f_shock(self.WL,p)
		else: return self.f_rarefaction(self.WL,p)
		
	def fR(self,p):
		if p>self.WR.p and self.WR.p != 0.: return self.f_shock(self.WR,p)
		else: return self.f_rarefaction(self.WR,p)
		
	def pressureFunc(self,p): #eqn (4.5) in Toro
		return self.fL(p) + self.fR(p) + (self.WR.u - self.WL.u)
		
	### Set p_star and u_star ###					
	def solve(self):
		#Numerically solve pressureFunc for p_star using Brent's method; get u_star from p_star
		self.p_star = brentq(self.pressureFunc,1e-40,1e20,rtol=1e-12,maxiter=200,full_output=False,disp=True)
		self.u_star = 0.5*(self.WL.u + self.WR.u) \
						+ 0.5*(self.fR(self.p_star)-self.fL(self.p_star))
	
	### Determine the left and right waves ###
	def waveStructure(self):
		if self.p_star != None:
			if self.p_star > self.WL.p:
				self.leftWave = 'shock'
			else:
				self.leftWave = 'rarefaction'
			if self.p_star > self.WR.p:
				self.rightWave = 'shock'
			else:
				self.rightWave = 'rarefaction'
		#Check for vacuum states
		if self.WL.p == 0.: 
			self.leftWave  = 'vacuum front'
			self.rightWave = 'rarefaction'
		if self.WR.p == 0.: 
			self.rightWave = 'vacuum front'
			self.leftWave = 'rarefaction'
		if self.p_star == 0.: #generated vacuum
			self.rightWave = 'rarefaction'
			self.leftWave = 'rarefaction'
		return {'leftWave':self.leftWave, 'rightWave':self.rightWave}
	
	### Solve density across the contact discontinuity ###
	def contactDisc(self,W,wave):
		if wave == 'shock':
			gfac =  (W.g-1.)/(W.g+1.)
			top = self.p_star/W.p + gfac
			bot = gfac*self.p_star/W.p + 1.
			return W.d*top/bot
		elif wave == 'rarefaction':
			return W.d*(self.p_star/W.p)**(1./W.g)
	
	### Gather variables in the star region ###
	def starVariables(self):
		self.dL_star = self.contactDisc(self.WL,self.leftWave)
		self.dR_star = self.contactDisc(self.WR,self.rightWave)
		starvals = {'dL_star':self.dL_star,'dR_star':self.dR_star,\
					'u_star':self.u_star,'p_star':self.p_star}
		return starvals
					
	### Tie everything together; print and record the solution ###
	def solution(self,display_output=True):
		soln = []
		ambientStates = {'left':self.WL,'right':self.WR}
		if (self.vacuum()): #no star region if vacuum region
			if self.WL.p == 0. or self.WR.p == 0.: self.p_star = None
			else: self.p_star = 0. #indicates vacuum is generated
			starRegion = {'dL_star':None,'dR_star':None,\
					'u_star':None,'p_star':self.p_star}
			waveStructure = self.waveStructure()
			if display_output:
				print "\n*** Solution to Riemann Problem: ***"
				print "Star Region: Does not exist in the presence of vacuum."
		else:
			self.solve() #solve the Riemann problem to get p_star,u_star
			waveStructure = self.waveStructure() #need wave structure to find density jump
			starRegion = self.starVariables() #determine density jump across contact
			if display_output:
				print "\n*** Solution to Riemann Problem: ***"
				print "Star Region:"
				if self.WL.g == 1. and self.WR.g == 1.:
					print "density = {}".format(starRegion['dL_star'])
				else:
					print "density left of contact = {}".format(starRegion['dL_star'])
					print "density right of contact = {}".format(starRegion['dR_star'])
				print "velocity = {}".format(starRegion['u_star'])
				print "pressure = {}".format(starRegion['p_star'])
		soln.append(ambientStates)
		soln.append(waveStructure) #record wave structure
		soln.append(starRegion) #calc density jump and save star variables
		if display_output:
			print "\nWave Structure: "
			print "left wave = {}".format(waveStructure['leftWave'])
			print "right wave = {}".format(waveStructure['rightWave'])
		return soln

# Class to separately hold, set, and then combine left and right state data profiles
class StateVector:
	def __init__(self,xL,xR):
		self.L = len(xL)
		self.R = len(xR)
		
	def set_scalarLstate(self,W=None):
		if W == None:
			self.gL = []
			self.dL = []
			self.uL = []
			self.pL = []
		else: 
			self.gL = np.ones(self.L)*W.g
			self.dL = np.ones(self.L)*W.d
			self.uL = np.ones(self.L)*W.u
			self.pL = np.ones(self.L)*W.p
		
	def set_scalarRstate(self,W=None):
		if W == None:
			self.gR = []
			self.dR = []
			self.uR = []
			self.pR = []
		else:
			self.gR = np.ones(self.R)*W.g
			self.dR = np.ones(self.R)*W.d
			self.uR = np.ones(self.R)*W.u
			self.pR = np.ones(self.R)*W.p

	def set_vectorLstate(self,W):
		self.gL = W[0]
		self.dL = W[1]
		self.uL = W[2]
		self.pL = W[3]
		
	def set_vectorRstate(self,W):
		self.gR = W[0]
		self.dR = W[1]
		self.uR = W[2]
		self.pR = W[3]
	
	def combineStates(self):
		g = np.concatenate((self.gL,self.gR))
		d = np.concatenate((self.dL,self.dR))
		u = np.concatenate((self.uL,self.uR))
		p = np.concatenate((self.pL,self.pR))
		W = []
		W.append(g); W.append(d); W.append(u); W.append(p);
		return W
		
# Class to sample the solution to the Riemann problem at a given (x,t) so as to 
# construct position profiles of an input variable at any given instant.  
class Exact1DProfile:
	def __init__(self,Solution,xaxis,time,center=0.,resize=True):
		ambientStates = Solution[0]
		waveStructure = Solution[1]
		starVariables = Solution[2]
		self.WL = ambientStates['left']
		self.WR = ambientStates['right']
		self.leftWave = waveStructure['leftWave']
		self.rightWave = waveStructure['rightWave']
		self.dL_star = starVariables['dL_star']
		self.dR_star = starVariables['dR_star']
		self.u_star = starVariables['u_star']
		self.p_star = starVariables['p_star']
		self.xaxis = xaxis
		self.center = center 
		if self.center != 0.: self.xaxis -= self.center
		if time==0. and (np.min(self.xaxis) >= 0. or np.max(self.xaxis) <=0.):
			print "WARNING: xaxis does not contain initial interface between ambient states!"
		self.t = time
		self.resize = resize
	
	### Slice axis on either side of xmark = wave_speed*time ###
	#If self.resize is True (default), the wavefront location, xmark, is added to both
	#split arrays; this is needed for the profiles to possess true discontinuities.
	#Preventing this resizing just means the exact wavefront won't appear to be tracked in
	#the resulting profile; this may desired if xaxis is fed in from a simulation, e.g.,
	#allowing an exact comparison with equally sized arrays.
	#(Note that xmark is 0 at time=0 and this 0 should not be appended.)
	def splitAxis(self,xaxis,speed):
		xL = []; xR = [];
		xmark = speed*self.t 
		if self.resize and self.t>0.: xR.append(xmark) 
		for x in xaxis:
			if x <= xmark:
				xL.append(x)
			else:
				xR.append(x)
		if self.resize and self.t>0.: xL.append(xmark)
		return np.array(xL),np.array(xR)

	### Combine separated left and right states ###
	def combineStates(self,WLvec,WRvec):
		g = np.concatenate((WLvec[0],WRvec[0]))
		d = np.concatenate((WLvec[1],WRvec[1]))
		u = np.concatenate((WLvec[2],WRvec[2]))
		p = np.concatenate((WLvec[3],WRvec[3]))
		W = []
		W.append(g); W.append(d); W.append(u); W.append(p);
		return W

	### Define relevant wave speeds ###
	def shockSpeed(self,W,state):
		d = W.d; u = W.u; p = W.p; a = W.a;
		gp1 = W.g+1.
		gm1 = W.g-1.
		ginv = 1./W.g
		if state == 'left':
			a *= -1.
		bracket = 0.5*gp1*ginv*self.p_star/p + 0.5*gm1*ginv
		return u + a * np.sqrt(bracket)
		
	def headSpeed(self,W,state):
		d = W.d; u = W.u; p = W.p; a = W.a;
		if state == 'left':
			a *= -1.
		return u + a

	def tailSpeed(self,W,state):
		d = W.d; u = W.u; p = W.p; a = W.a;
		gm1 = W.g-1.
		a_star = a*(self.p_star/p)**(0.5*gm1/W.g)
		if state == 'left':
			a_star *= -1.
		return self.u_star + a_star
		
	def vacuumFrontSpeed(self,W,state):
		d = W.d; u = W.u; p = W.p; a = W.a;
		if W.g == 1.:  #ISOTHERMAL CASE
			hfac = np.log(d)
		else: hfac = 2./(W.g-1.)
		if state == 'left':
			a *= -1.
		return u - hfac*a
		
	### Find self-similar solution in a rarefaction fan ###
	def fanRegion(self,W,x,state):
		d = W.d; u = W.u; p = W.p; a = W.a;
		ones = np.ones(len(x))
		t = self.t
		if t==0.:
			t+=1e-10
		if state == 'left':
			a *= -1.
		gm1 = W.g-1.
		gp1_inv = 1./(W.g+1.)
		ufan = 2.*gp1_inv*(-a*ones + 0.5*gm1*u*ones + x/t)
		if W.g == 1.0: #ISOTHERMAL CASE
			dfan = d*np.exp(-u*ones/a + x/(a*t) - ones)
			pfan = a*a*dfan
		else: #ADIABATIC CASE
			dfan = 2.*gp1_inv*ones - gm1*gp1_inv*(u*ones-x/t)/a
			#clear any zeros which can occur if resize == True in vacuum cases
			dfan = [val if val>0. else 0. for val in dfan]
			#ignore tiny numbers which occur if resize == True in vacuum cases
			dfan = [d*val**(2./gm1) if np.abs(val) > 1e-12 else val for val in dfan]
			pfan = p*(2.*gp1_inv*ones - gm1*gp1_inv*(u*ones-x/t)/a)**(2.*W.g/gm1)
		Wfan = []
		Wfan.append(W.g*ones)
		Wfan.append(dfan)
		Wfan.append(ufan)
		Wfan.append(pfan)
		return Wfan
	
	### Construct profiles if left/right wave is a shock ###
	def shock(self,x,state):
		if state == 'left':
			v_shock = self.shockSpeed(self.WL,state)
			xL,xR = self.splitAxis(x,v_shock)
			Wvec = StateVector(xL,xR)
			Wvec.set_scalarLstate(self.WL)
			Wstar = State(self.WL.g,self.dL_star,self.u_star,self.p_star)
			Wvec.set_scalarRstate(Wstar)
		else:
			v_shock = self.shockSpeed(self.WR,state)
			xL,xR = self.splitAxis(x,v_shock)
			Wvec = StateVector(xL,xR)
			Wvec.set_scalarRstate(self.WR)
			Wstar = State(self.WR.g,self.dR_star,self.u_star,self.p_star)
			Wvec.set_scalarLstate(Wstar)
		xnew = np.concatenate((xL,xR))
		return xnew,Wvec.combineStates()
		
	### Construct profiles if left/right wave is a rarefaction ###
	def rarefaction(self,x,state):
		if state == 'left':
			v_head = self.headSpeed(self.WL,state)
			xL,xblank = self.splitAxis(x,v_head)
			Wvec1 = StateVector(xL,[])
			Wvec1.set_scalarLstate(self.WL)
			Wvec1.set_scalarRstate()
			W1 = Wvec1.combineStates()
			v_tail = self.tailSpeed(self.WL,state)
			xbounds = np.array([v_head*self.t,v_tail*self.t])
			Wbounds = self.fanRegion(self.WL,xbounds,state)
			xfan,xR = self.splitAxis(xblank,v_tail)
			Wvec2 = StateVector(xfan,xR)
			Wfan = self.fanRegion(self.WL,xfan,state)
			Wvec2.set_vectorLstate(Wfan)
			Wstar = State(self.WL.g,self.dL_star,self.u_star,self.p_star)
			Wvec2.set_scalarRstate(Wstar)
			W2 = Wvec2.combineStates()
			W = self.combineStates(W1,W2)
		else:
			v_head = self.headSpeed(self.WR,state)
			xblank,xR = self.splitAxis(x,v_head)
			Wvec1 = StateVector([],xR)
			Wvec1.set_scalarRstate(self.WR)
			Wvec1.set_scalarLstate()
			W1 = Wvec1.combineStates()
			v_tail = self.tailSpeed(self.WR,state)
			xbounds = np.array([v_tail*self.t,v_head*self.t])
			Wbounds = self.fanRegion(self.WR,xbounds,state)
			xL,xfan = self.splitAxis(xblank,v_tail)
			Wvec2 = StateVector(xL,xfan)
			Wfan = self.fanRegion(self.WR,xfan,state)
			Wvec2.set_vectorRstate(Wfan)
			Wstar = State(self.WR.g,self.dR_star,self.u_star,self.p_star)
			Wvec2.set_scalarLstate(Wstar)
			W2 = Wvec2.combineStates()
			W = self.combineStates(W2,W1)
		xnew = np.concatenate((xL,xfan,xR))
		return xnew,W
				
	### Construct profiles in the case of vacuum ###
	def vacuum(self,x,state):
		if state == 'left': #left rarefaction if RIGHT vacuum region
			xL,xblank = self.splitAxis(x,self.WL.u - self.WL.a)
			Wvec1 = StateVector(xL,[])
			Wvec1.set_scalarLstate(self.WL)
			Wvec1.set_scalarRstate()
			W1 = Wvec1.combineStates()
			v_front = self.vacuumFrontSpeed(self.WL,state)
			xfan,xR = self.splitAxis(xblank,v_front)
			Wvec2 = StateVector(xfan,xR)
			Wfan = self.fanRegion(self.WL,xfan,state)
			Wvec2.set_vectorLstate(Wfan)
			Wvacuum = State(self.WL.g,0.,v_front,0.)
			Wvec2.set_scalarRstate(Wvacuum)
			W2 = Wvec2.combineStates()
			W = self.combineStates(W1,W2)
			xnew = np.concatenate((xL,xfan,xR))
		elif state == 'right': #right rarefaction if LEFT vacuum region
			xblank,xR = self.splitAxis(x,self.WR.u + self.WR.a)
			Wvec1 = StateVector([],xR)
			Wvec1.set_scalarRstate(self.WR)
			Wvec1.set_scalarLstate()
			W1 = Wvec1.combineStates()
			v_front = self.vacuumFrontSpeed(self.WR,state)
			xL,xfan = self.splitAxis(xblank,v_front)
			Wvec2 = StateVector(xL,xfan)
			Wfan = self.fanRegion(self.WR,xfan,state)
			Wvec2.set_vectorRstate(Wfan)
			Wvacuum = State(self.WR.g,0.,v_front,0.)
			Wvec2.set_scalarLstate(Wvacuum)
			W2 = Wvec2.combineStates()
			W = self.combineStates(W2,W1)
			xnew = np.concatenate((xL,xfan,xR))
		elif state == 'center': #generated vacuum case
			v_frontR = self.vacuumFrontSpeed(self.WR,'right')
			xblankL,xblankR = self.splitAxis(x,v_frontR) #
			xL,xblankL2 = self.splitAxis(xblankL,self.WL.u - self.WL.a)
			Wvec1 = StateVector(xL,[])
			Wvec1.set_scalarLstate(self.WL)
			Wvec1.set_scalarRstate()
			W1 = Wvec1.combineStates()
			
			v_frontL = self.vacuumFrontSpeed(self.WL,'left')
			xfanL,xvacuum = self.splitAxis(xblankL2,v_frontL)
			Wvec2 = StateVector(xfanL,xvacuum)
			WfanL = self.fanRegion(self.WL,xfanL,'left')
			
			Wvec2.set_vectorLstate(WfanL)
			Wvacuum = State(self.WL.g,0.,v_frontL,0.) #CHECK v_front
			Wvec2.set_scalarRstate(Wvacuum)
			W2 = Wvec2.combineStates()
			WL = self.combineStates(W1,W2)
			
			xfanR,xR = self.splitAxis(xblankR,self.WR.u + self.WR.a)
			Wvec3 = StateVector(xfanR,xR)
			WfanR = self.fanRegion(self.WR,xfanR,'right')
			Wvec3.set_vectorLstate(WfanR)
			Wvec3.set_scalarRstate(self.WR)
			WR = Wvec3.combineStates()
			W = self.combineStates(WL,WR)
			xnew = np.concatenate((xL,xfanL,xvacuum,xfanR,xR))
		return xnew,W
		
	### Build profiles based on solution to Riemann problem ###
	def Wprofiles(self):
		if self.p_star == None: #L/R vacuum state case
			if self.WR.p == 0.:
				xnew,W = self.vacuum(self.xaxis,'left')
			elif self.WL.p == 0.: 
				xnew,W = self.vacuum(self.xaxis,'right')
		elif self.p_star == 0.: #generated vacuum case
			xnew,W = self.vacuum(self.xaxis,'center')
		else: #non-vacuum case
			xL,xR = self.splitAxis(self.xaxis,self.u_star)
			if self.leftWave == 'shock':
				xLnew,WLvec = self.shock(xL,'left')
			else: 
				xLnew,WLvec = self.rarefaction(xL,'left')
			if self.rightWave == 'shock':
				xRnew,WRvec = self.shock(xR,'right')
			else: 
				xRnew,WRvec = self.rarefaction(xR,'right')
			xnew = np.concatenate((xLnew,xRnew))
			W = self.combineStates(WLvec,WRvec)
		return xnew,W
	
	def getVariables(self,choice=None,tails=False,custom=False):
		vars = []
		if custom: #then choice must be a list of strings
			for var in choice:
				vars.append(var)
		else:
			if choice == None: #then all variables are output
				vars.append('x')
				vars.append('density')
				vars.append('velocity')
				vars.append('pressure')
				vars.append('momentum')
				vars.append('sound speed')
				vars.append('Mach number')
				vars.append('internal energy')
				vars.append('total energy')
			elif choice == 'primitive': #then primitive variables are output
				vars.append('x')
				vars.append('density')
				vars.append('velocity')
				vars.append('pressure')
			elif choice == 'conserved': #then conserved variables are output
				vars.append('x')
				vars.append('density')
				vars.append('momentum')
				vars.append('total energy')
			elif choice == 'other': #then remaining variables are output
				vars.append('x')
				vars.append('sound speed')
				vars.append('Mach number')
				vars.append('internal energy')
			elif choice[0] == 'comoving': #comoving frame data is output
				vars.append(choice[1])
				if self.u_star == None:
					vars.append('x')
				if self.leftWave == 'shock':
					vars.append('(left) shock frame')
				elif self.leftWave == 'rarefaction': 
					if tails:
						vars.append('(left) rarefaction tail frame')
					else:
						vars.append('(left) rarefaction head frame')
				else:
					vars.append('(left) vacuum front frame')
				if self.rightWave == 'shock':
					vars.append('(right) shock frame')
				elif self.rightWave == 'rarefaction': 
					if tails:
						vars.append('(right) rarefaction tail frame')
					else:
						vars.append('(right) rarefaction head frame')
				else:
					vars.append('(right) vacuum front frame')
				if self.u_star is not None:
					vars.append('contact discontinuity frame')
			else:
				print "\nERROR: choice \"{}\" not recognized!\n".format(choice)
				print "Choose from defaults:"\
						"\nChoice = 'primitive': plots density, velocity, and pressure"\
						"\nChoice = 'conserved': plots density, momentum, and total energy"\
						"\nChoice = 'other': plots sound speed, Mach number, and internal energy"\
						"\nChoice = ['comoving', var]: plots var in the comoving frame of left/center/right waves"\
						"\nFor the comoving choices, var should be chosen from:\n"\
						"'x' 'density', 'velocity', 'pressure', 'sound speed',\n" \
						"'momentum', 'Mach number', 'internal energy', 'total energy'\n" \
						"\nCalling run with run(choice,tail=True) uses the frame of the tail of any rarefaction waves instead of the head.\n"\
						"\nAlternatively, is custom is set to True, set choice to a list of "\
						"strings with any of the above (dependent) variables.\n"\
						"You can also choose different comoving frame "\
						"(independent) variables from: \n"\
						"'(left) rarefaction head frame', '(right) rarefaction head frame'\n"\
						"'(left) rarefaction tail frame', '(right) rarefaction tail frame'\n"\
						"'(left) shock frame', '(right) shock frame'\n"\
						"'(left) vacuum front frame', '(right) vacuum front frame'\n"\
						"'contact discontinuity frame'\n"
				sys.exit()
		return vars
		
	### Tie everything together; build default or custom output profiles ###
	def run(self,choice=None,tails=False,custom=False):
		x,W = self.Wprofiles()
		if self.center != 0.: x += self.center #return x to its original interval
		gamma = W[0]
		density = W[1]		
		velocity = W[2]			
		pressure = W[3]
		blank = np.ones(len(x))
		#change zeros to nans in density
		density = np.array([val if np.abs(val) > 1e-16 else np.nan for val in density])
		
		profiles = {}
		vars = self.getVariables(choice,custom)
		for var in vars:
			if var == 'x':
				y = x
			elif var == 'density':
				y = density
			elif var == 'velocity':
				y = velocity
			elif var == 'pressure':
				y = pressure 
			elif var == 'momentum':
				y = density*velocity
			elif var == 'sound speed':
				y = np.sqrt(gamma*pressure/density)
			elif var == 'Mach number':
				y = velocity/np.sqrt(gamma*pressure/density)
			elif var == 'internal energy':
				if self.WL.g == 1. or self.WR.g == 1.:
					gamma += 1e-14
				y = pressure/(density*(gamma-1.))	
			elif var == 'total energy':
				if self.WL.g == 1. or self.WR.g == 1.:
					gamma += 1e-14
				y = pressure/(density*(gamma-1.)) + 0.5*density*velocity**2
			elif var == '(left) rarefaction head frame':
				v_head = self.headSpeed(self.WL,'left')
				y = x - v_head*blank*self.t
			elif var == '(left) rarefaction tail frame':
				if self.u_star == None: #vacuum cases
					v_tail = self.vacuumFrontSpeed(self.WL,'left')
				else:
					v_tail = self.tailSpeed(self.WL,'left')
				y = x - v_tail*blank*self.t
			elif var == '(left) shock frame':
				gp1 = self.WL.g+1.; gm1 = self.WL.g-1.; ginv = 1./self.WL.g;
				a = -np.sqrt(self.WL.g*self.WL.p/self.WL.d)
				bracket = np.sqrt(0.5*gp1*ginv*self.p_star/self.WL.p + 0.5*gm1*ginv)
				y = x - (self.WL.u + a*bracket)*blank*self.t
			elif var == '(left) vacuum front frame':
				y = x - self.vacuumFrontSpeed(self.WR,'right')*blank*self.t
			elif var == 'contact discontinuity frame':
				y = x - self.u_star*blank*self.t			
			elif var == '(right) rarefaction head frame':
				v_head = self.headSpeed(self.WR,'right')
				y = x - v_head*blank*self.t
			elif var == '(right) rarefaction tail frame':
				if self.u_star == None: #vacuum cases
					v_tail = self.vacuumFrontSpeed(self.WR,'right')
				else:
					v_tail = self.tailSpeed(self.WR,'right')
				y = x - v_tail*blank*self.t
			elif var == '(right) shock frame':
				gp1 = self.WR.g+1.; gm1 = self.WR.g-1.; ginv = 1./self.WR.g;
				a = np.sqrt(self.WR.g*self.WR.p/self.WR.d)
				bracket = np.sqrt(0.5*gp1*ginv*self.p_star/self.WR.p + 0.5*gm1*ginv)
				y = x - (self.WR.u + a*bracket)*blank*self.t
			elif var == '(right) vacuum front frame':
				y = x - self.vacuumFrontSpeed(self.WL,'left')*blank*self.t
			else: 
				print "\nERROR: variable \"{}\" not recognized!\n".format(var)
				print "Choose from:\n"\
						"'x' 'density', 'velocity', 'pressure', 'sound speed'," \
						"'momentum', 'Mach number', 'internal energy', 'total energy'\n" \
						"'(left) rarefaction head frame', '(left) rarefaction tail frame'\n"\
						"'(right) rarefaction head frame', '(right) rarefaction tail frame'\n"\
						"'(left) shock frame', '(right) shock frame'\n"\
						"'(left) vacuum front frame', '(right) vacuum front frame'\n"\
						"'contact discontinuity frame'\n"
			profiles[var] = y
		return profiles
		
	### Determine max/min of variables so can auto set plot limits ###
	def getLimits(self,pad=1.1):
		profiles = self.run()
		
		#Clean out any nans or infs which occur for vacuum cases
		d = [x for x in profiles['density'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		a = [x for x in profiles['sound speed'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		m = [x for x in profiles['momentum'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		M = [x for x in profiles['Mach number'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		Egas = [x for x in profiles['internal energy'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		Etot = [x for x in profiles['total energy'] if np.abs(x)!=np.inf and np.isnan(x) == False]
		
		#Energy for isothermal vacuum cases is all nans and infs leaving empty arrays. 
		if len(Egas) == 0:
			Egas = [0]
		if len(Etot) == 0:
			Etot = [0]
		
		mins = {'density':	min(0.,np.min(d)),\
				'velocity':	min(0.,pad*np.min(profiles['velocity'])),\
				'pressure':	min(0.,np.min(profiles['pressure'])),\
				'momentum':	min(0.,pad*np.min(m)),\
				'sound speed':	min(0.,np.min(a)),\
				'Mach number':	min(0.,pad*np.min(M)),\
				'internal energy':	min(0.,np.min(Egas)),\
				'total energy':	min(0.,np.min(Etot))}
		maxs = {'density':	pad*np.max(d),\
				'velocity':	1.+pad*np.max(profiles['velocity']),\
				'pressure':	pad*np.max(profiles['pressure']),\
				'momentum':	pad*np.max(m),\
				'sound speed':	pad*np.max(a),\
				'Mach number':	pad*np.max(M),\
				'internal energy':	pad*np.max(Egas),\
				'total energy':	pad*np.max(Etot)}
		
		#Offset max by pad if max is zero		
		for var in maxs:
			if maxs[var] == 0.:
				maxs[var] += pad
		
		limits = []
		limits.append(mins)
		limits.append(maxs)
		return limits

"""
--------- EXAMPLE SCRIPT ---------
Uncomment the section below to plot the solution to the classic Sod shock-tube test.
By default, inputting one variable will output position profiles of that variable as a
function of x (considered the lab frame - top left plot) as well as 3 comoving frames:
1. The frame of the contact discontinuity (top right plot)
2. The frame of the left rarefaction wave (bottom left plot)
3. The frame of the right shock wave (bottom right plot)  
The wavefronts of each of these waves is stationary in their comoving frames. Therefore,
overplotting simulation data would allow one to test whether the hydrodynamical scheme
is Galilean invariant, as it should be.  

Note that the left and right wave-structure is not known a priori; a shock occurs if the
pressure in the unknown star region exceeds the ambient pressure.  Otherwise the wave
is a rarefaction.  Also note that the isothermal (gamma = 1.0) Riemann problem does not
have a contact discontinuity.  (The governing equations have two eigenvectors rather 
than three; the third arises from the energy equation, which is not present in the 
isothermal case as energy is no longer conserved.)

By default, inputting one variable will output position profiles of that variable as a
function of x (considered the lab frame - top left plot) as well as 3 comoving frames:
1. The frame of the contact discontinuity (top right plot)
2. The frame of the left rarefaction wave (bottom left plot)
3. The frame of the right shock wave (bottom right plot)  
The wavefronts of each of these waves is stationary in their comoving frames. Therefore,
overplotting simulation data would allow one to test whether the hydrodynamical scheme
is Galilean invariant, as it should be.

To run other Riemann problems, simply change the input states WL and WR. 
"""

"""
#---------------------------------------------------------------------
WL = State(1.4,1.0,0.0,1.0)   #gamma=1.4, density=1.0,   velocity=0., pressure=1.0
WR = State(1.4,0.125,0.0,0.1) #gamma=1.4, density=0.125, velocity=0., pressure=0.1

xmin =-0.5
xmax = 0.5
xc = 0.
Nx = 1000
xargs = [xmin,xmax,Nx,xc]
xaxis = np.linspace(xmin,xmax,Nx)

#choice = 'primitive'
#choice = 'conserved'
#choice = 'other'
choice = ['comoving','velocity']

N = 200
delay = 100
DEBUG_MODE = False

#-----------------------------------------------------------------
problem = RiemannProblem(WL,WR)
solution = problem.solution()

xmin = np.min(xaxis); xmax = np.max(xaxis);
query = Exact1DProfile(solution,xaxis,0.01,center=0.)
limits = query.getLimits()
min = limits[0]; max = limits[1] #dictionaries of plot limits

profiles = query.run(choice)
y = query.getVariables(choice) #list of 3 variables

if DEBUG_MODE:
	print "\n\nMinimums = {} \n\nMaximums = {}".format(min,max)
	for var in profiles:
		print "\n\n{} = ".format(var),profiles[var]

#crossingTime = 0.5*(xmax-xmin)/abs(WR.u-WL.u + 0.5*(WL.a+WR.a))
crossingTime = (xmax-xmin)/np.max([abs(WR.u), WR.a, abs(WL.u), WL.a])
dt = crossingTime/N
fig = plt.figure(figsize=(18, 6))
L = fig.add_subplot(131); M = fig.add_subplot(132); R = fig.add_subplot(133); 
	
L.set_xlim(xmin, xmax); 
Lline, = L.plot([], [], 'r', lw=2)
M.set_xlim(xmin, xmax); 
Mline, = M.plot([], [], 'r', lw=2)
R.set_xlim(xmin, xmax); 
Rline, = R.plot([], [], 'r', lw=2)

if choice[0] == 'comoving':
	L.set_ylim(min[y[0]], max[y[0]]);
	L.set_xlabel(y[1]); L.set_ylabel(y[0])	
	M.set_ylim(min[y[0]], max[y[0]]);
	M.set_xlabel(y[3]);
	R.set_ylim(min[y[0]], max[y[0]]);
	R.set_xlabel(y[2]);
else:
	L.set_ylim(min[y[1]], max[y[1]]);
	L.set_xlabel('x'); L.set_ylabel(y[1])	
	M.set_ylim(min[y[2]], max[y[2]]);
	M.set_xlabel('x'); M.set_ylabel(y[2])	
	R.set_ylim(min[y[3]], max[y[3]]);
	R.set_xlabel('x'); R.set_ylabel(y[3])

#Display time above TL window
time_template = 'time = %.2f'
time_text = L.text(0.05, 0.94, '', transform=L.transAxes)

def demo_plot(time,xmin,xmax,Nx,xc):
	
	xaxis = np.linspace(xmin,xmax,Nx)
	analyticsoln = Exact1DProfile(solution,xaxis,time,center=xc)
	profiles = analyticsoln.run(choice)
	
	# Load profiles
	x = profiles[y[0]]
	y1 = profiles[y[1]]
	y2 = profiles[y[2]]
	y3 = profiles[y[3]]	
	
	# Plot profiles
	if choice[0] == 'comoving':
		Lline.set_data(y1,x);
		Rline.set_data(y2,x);
		Mline.set_data(y3,x);
	else:
		Lline.set_data(x,y1);
		Mline.set_data(x,y2);
		Rline.set_data(x,y3);	
	# Annotate time on screen
	time_text.set_text(time_template%(time))

time = 0.1
demo_plot(time,xmin,xmax,Nx,xc)
plt.show()
"""