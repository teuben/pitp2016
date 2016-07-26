The following packages we are likely to need at PITP 2016. Each UNIX
distro (linux and mac alike) have different ways how to install
them. Some example are given here how to do them.

You can find copies of codes on either our published URL's , or
within the directory /data/pitp on the IAS cluster where all participants
have an account. Look in /data/pitp/software.
If you want to work on one of the cluster nodes, be sure to create your
own directory within /data/pitp/ and work within that.
Please read http://www.sns.ias.edu/computing/servers carefully where you
should and should not compute. Use the commands qsub/qstat/qalter/qresub
to manage jobs that need resources.


Example copying files from laptop to cluster
    scp readme pitp80@ssh.sns.ias.edu:/data/pitp/pitp80

Example copying file from cluster to laptop :
    scp pitp80@ssh.sns.ias.edu:/data/pitp/pitp80/run1.dat .

--------------------------------------------------------------------------------
UNIX COMMANDS:


## git

  we also use this to distribute notes, scripts, software for PITP 2016.  
  To get a copy of the PITP notes:

    git clone http://github.com/teuben/pitp2016
    cd pitp2016

  and to get a updated version with the latest additions:

    cd pitp2016
    git pull

  See also Lupton's talk and his github repo on:
  https://github.com/RobertLuptonTheGood/APC524GitLecture


## gcc/g++/gfortran

  The C/C++/Fortran(90) compiler suite.  Others may be available as well that work
  also, e.g. "clang", the Intel compiler and the PGI compiler.  
  The C compiler needs to suport OpenMP (which gcc does)


## python/ipython

  Although the python command may be present for you, what we really want here is
  an integrated scientific python environment, e.g. via anaconda.  
  Within python, the following modules are recommended:

	numpy
	scipy
	matplotlib
	yt
	astropy
	pandas
	scikit-learn
	glueviz


## ffmpeg

  This tool creates movies from still images that simulations often create on the fly.  
  Great for presentations. We will use this for a few of the homeworks.

  In ubuntu linux:
    
    sudo apt install ffmpeg imagemagick vlc

## gnuplot

  A simple, basic but effective way to plot data.

## ygraph

  Another simple and basic graph plotter.   
  See for example http://cactuscode.org/documentation/visualization/yGraph/

## visit

  A more graphical based visualization environment  
  See https://wci.llnl.gov/simulation/computer-codes/visit/


--------------------------------------------------------------------------------
UNIX LIBRARIES:

## openmpi 

  This is one of the ways to implement MPI, and probably our preferred way.  
  On Ubuntu/Linux:
  
    apt install openmpi-bin libopenmpi-dev libhdf5-openmpi-dev

  On mac/homebrew:
    
    brew install gcc
    brew install openmpi --enable-mpi-thread-multiple
    brew install hdf5 --with-fortran --with-mpi

## mpich

  This is another MPI library:  
  On Ubuntu/Linux:
  
    apt install mpich
  
  but not experimented with this yet, as openmpi seems to work.

## gsl

  This is the Gnu Scientific Library. Use your package manager to install it, e.g.  
  In Ubuntu linux:
            
    sudo apt install libgsl-dev

## yt  (ZuHone)
	
  Versatile visualization environment in python. Probably easiest to install
  within your python environment via the commands
  
    conda install yt
    or 
    pip install yt
	
  depending which type of python environment you have.
  
  Lecture notes are here:  http://yt-project.org/pitp2016_demo/yt_tutorial.slides.html

  The new yt 3.3-dev (development) series can be installed as follows:  

    conda install -c http://use.yt/with_conda yt                 (unix)
    conda install -c http://conda.anaconda.org/jzuhone yt        (windows)
		
  or if you have another python stack (or just want to use pip):
  
    # build from source on macOS or Linux:
    pip install hg+https://bitbucket.org/yt_analysis/yt#egg=yt
    # binary install on macOS or Windows:
    pip install --pre yt

	wget -r --no-parent http://yt-project.org/pitp2016_demo/data/WindTunnel
	wget -r --no-parent http://yt-project.org/pitp2016_demo/data/DD0046
	wget -r --no-parent http://yt-project.org/pitp2016_demo/examples


## d3js  (Ericson)

  This is a JavaScript based visualization toolkit. Great to present your research online
  for outreach, collaborators etc.  
  See https://d3js.org/


## athena  (Stone)

  http://www.astro.princeton.edu/~jstone/downloads/papers/stone_hw.pdf  
  Use athena4.2.tar.gz

## harmpi  (Tchekhovskoy)

  GitHub repo: https://github.com/atchekho/harmpi.git  
  Tutorial:    https://github.com/atchekho/harmpi/blob/master/tutorial.md  
  Homework:    https://github.com/atchekho/harmpi/blob/master/exercises.md  

  The code is 3D and is parallelized using MPI, results can be obtained rather quickly.

## iharm2d_v3 (Gammie)

  GitHub repo 1: https://github.com/AFD-Illinois/iharm2d_v3

  GitHub repo 2: https://github.com/AFD-Illinois/ibothros2d

  see also additional files in pitp2016/Gammie

  For some version of mac (e.g. using brew) you may find that gsl has been installed in /usr/local
  instead of /usr.  In that case the -I/usr/local/include and -L/usr/local/lib flags may need to
  be added to the Makefile CFLAGS and LDFLAGS macros resp.
  Also, for some of those I've seen the compiler is called gcc-6, so change CC=gcc-6 to fix the
  the compilation

## tristan (Spitovsky) 

  Homework at:   http://www.astro.princeton.edu/~anatoly/PICHW2016.html  
  Code at:       https://github.com/ntoles/tristan-mp-pitp  
  Vis at:        https://github.com/pcrumley/Iseult

    <edit> Makefile
    make
    mpirun -n 1 ./tristan-mp2d -i input.weibel
    Iseult/iseult.py

  Here it is important to have an MPI enabled compiler !!! See the "openmpi"
  comments earlier in this document.

## ADERG (Zanotti)

   See code in pitp2016/Zanotti
