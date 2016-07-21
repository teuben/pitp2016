"""
--------- DESCRIPTION ---------
This script is an example of how to animate 1D simulation
output from Athena4.2 (tab files).  There are two external
programs imported:
1. AthenaTabData.py - the actual reader for 1D Athena data
2. RiemannNRHD.py - a program to solve the Riemann problem
This script assumes you have run two different Athena
simulations using the shkset1D.c problem generator, one
at resolution 64 and another at resolution 256.  It then
animates these two solutions, along with the analytic 
solution, which is calculated by RiemannNRHD.py

If you encounter any issues, feel free to email me at 
waterst3@unlv.nevada.edu
"""

import array as ar
import numpy as np
import pylab as plt
import matplotlib
import AthenaTabData 
import RiemannNRHD as RP
from matplotlib import animation

# custom plotting routine
def plot_xy(line,x,y,c='k',m=None,ls=None,lw=1,label=''):
	
	# plot data
	line.set_data(x,y); 
	
	# use a solid line by default
	if m == None and ls == None: 
		ls = '-' 
	
	# set plot styles
	line.set_color(c)
	line.set_linewidth(lw)
	line.set_label(label)
	if m != None:
		line.set_marker(m); 
	if ls != None:
		line.set_linestyle(ls); 

#-----------------------------------------------------------------
""" input options """

ANIMATE = 1     # play animation on screen
MAKE_MOVIE = 1  # make a mp4 movie of animation

# pick simulation dumps to plot
dumps = ar.array('i',(i for i in range(0,81,1)))
Ndumps = len(dumps) #number of frames to animate

#-----------------------------------------------------------------
""" load athena simulations """
Sim1 = {}
Sim1['simpath'] = '' 
Sim1['simfolder'] = 'tst64'
Sim1['simtag'] = 'ToroSod'
Sim1['timestep'] = 0.005
Sim1['gamma'] = 1.4

Sim2 = {}
Sim2['simpath'] = '' 
Sim2['simfolder'] = 'tst256'
Sim2['simtag'] = 'ToroSod'
Sim2['timestep'] = 0.005
Sim2['gamma'] = 1.4

sim1 = AthenaTabData.AthenaSim1D(Sim1)
sim2 = AthenaTabData.AthenaSim1D(Sim2)

#-----------------------------------------------------------------
""" calculations """

# specify shock tube conditions
WL = RP.State(1.4,1.0,0.75,1.0)   #Left state: gamma, density, velocity, pressure
WR = RP.State(1.4,0.125,0.0,0.1)  #Right state: gamma, density, velocity, pressure

# spatial domain of analytic solution
xmin =-0.5; xmax = 0.5
xc = 0. # location of interface
Nx = 1000 # number of points 
xargs = [xmin,xmax,Nx,xc]
xaxis = np.linspace(xmin,xmax,Nx)

# choose which quantities to plot
choice = 'primitive'  # (density, velocity, pressure)
#choice = 'conserved' # (density, momentum, total energy)
#choice = 'other' # (sound speed, Mach number, and internal energy)

# solve the Riemann problem
problem = RP.RiemannProblem(WL,WR)
solution = problem.solution()

# query the solution to get y-limits
xmin = np.min(xaxis); xmax = np.max(xaxis);
query = RP.Exact1DProfile(solution,xaxis,0.01,center=0.)
limits = query.getLimits()
min = limits[0]; max = limits[1] #dictionaries of plot limits

profiles = query.run(choice)
y = query.getVariables(choice) #list of 3 variables

# set sim timestep and determine crossing time
dt = Sim1['timestep'] # note - assumes Sim2 has the same dt!
max_speed = np.max([abs(WR.u), WR.a, abs(WL.u), WL.a])
t_crossing = 0.5*(xmax-xmin)/max_speed
print "\ncrossing time = {}".format(t_crossing)

#-----------------------------------------------------------------
""" setup plot figure """

fig = plt.figure(figsize=(8, 8))

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = '14'
matplotlib.rcParams['ps.fonttype'] = 42 # compatible with MNRAS style file

A = fig.add_subplot(311)
B = fig.add_subplot(312)
C = fig.add_subplot(313)

# set plot limits	
A.set_xlim(xmin, xmax)
B.set_xlim(xmin, xmax)
C.set_xlim(xmin, xmax)
A.set_ylim(min[y[1]], max[y[1]])
B.set_ylim(min[y[2]], max[y[2]])
C.set_ylim(min[y[3]], max[y[3]])

# set plot labels
A.set_xlabel('x'); A.set_ylabel(y[1])
B.set_xlabel('x'); B.set_ylabel(y[2])
C.set_xlabel('x'); C.set_ylabel(y[3])

# title
A.set_title('Toro\'s Sod shock tube solution at different resolutions');

# initialize some plotting lines
a1,a2,a3 = A.plot([], [], [], [], [], [], 'r', lw=1)
b1,b2,b3 = B.plot([], [], [], [], [], [], 'r', lw=1)
c1,c2,c3 = C.plot([], [], [], [], [], [], 'r', lw=1)

#Display time above TL window
time_template = r'$t [t_{cross}] = %.3f$'
time_text = A.text(0.75, 0.85, '', transform=A.transAxes, size='large')

#--------------------------------------------------------------------- 
""" define animation function """

def animate(iter):
	time = dumps[iter]*dt
	
	# plot sim1
	sim1.load_data(dumps[iter])
	x_sim1 = sim1.hydro_profile('x')
	d_sim1 = sim1.hydro_profile('density')
	v_sim1 = sim1.hydro_profile('velocity')
	p_sim1 = sim1.hydro_profile('pressure')
	
	plot_xy(a1,x_sim1,d_sim1,c='r',lw=2,label='Nx = 64')
	plot_xy(b1,x_sim1,v_sim1,c='r',lw=2)
	plot_xy(c1,x_sim1,p_sim1,c='r',lw=2)
	
	# plot  sim2
	sim2.load_data(dumps[iter])
	x_sim2 = sim2.hydro_profile('x')
	d_sim2 = sim2.hydro_profile('density')
	v_sim2 = sim2.hydro_profile('velocity')
	p_sim2 = sim2.hydro_profile('pressure')
	
	plot_xy(a2,x_sim2,d_sim2,c='c',lw=2,label='Nx = 256')
	plot_xy(b2,x_sim2,v_sim2,c='c',lw=2)
	plot_xy(c2,x_sim2,p_sim2,c='c',lw=2)
	
	# plot analytic solution
	xaxis = np.linspace(xmin,xmax,Nx)
	analytic_soln = RP.Exact1DProfile(solution,xaxis,time,center=xc)
	profiles = analytic_soln.run(choice)
	x = profiles[y[0]]
	y1 = profiles[y[1]]
	y2 = profiles[y[2]]
	y3 = profiles[y[3]]	
	
	plot_xy(a3,x,y1,label='analytic')
	plot_xy(b3,x,y2)
	plot_xy(c3,x,y3)
		
	# Annotate time on screen
	time_text.set_text(time_template%(time/t_crossing))
	
	# Add legend
	legend = A.legend(loc='lower left', shadow=True, fontsize='large', fancybox=True, frameon=True) #, framealpha=0.6)

#--------------------------------------------------------------------- 
""" execute """

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=None, hspace=0.25)

if ANIMATE:
	ani = animation.FuncAnimation(fig,animate,init_func=None,
    frames=Ndumps,interval=100, blit=False,repeat=False)
else:
	animate(10)	
	
if ANIMATE and MAKE_MOVIE:
	#moviepath = '/Users/twaters/Movies/'
	name = 'ToroSod.mp4'
	#name = ''.join([i for i in name if i!='/']) #remove any '/' in name
	ani.save(name, fps=10, bitrate=1000) 
else:
	plt.show()

