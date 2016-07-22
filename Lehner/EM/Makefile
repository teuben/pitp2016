PROG =	wave

OBJS =	boundary.o derivs.o evolve.o initial.o main.o params.o rhs.o  

MAKEFLAGS = -r

LIBS = # -L/usr/local/lib/ -lbbhutil -lrnpl -lutilio -lvutil -lxvs 

F90 = ifort 
F77 = ifort 
F90FLAGS = -r8 -O3 -cpp 
# -g -fpstkchk -fpe0 -ftz
LDFLAGS = -L/opt/intel/fce/9.0/lib
LIBPATH = -L/usr/local/lib/

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.g90

.SUFFIXES: .o .mod .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f90.mod:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

.f.mod:
	$(F90) $(F90FLAGS) -c $<

boundary.o: params.o
derivs.o: params.o
evolve.o: boundary.o derivs.o params.o rhs.o
initial.o: derivs.o
main.o: derivs.o evolve.o initial.o params.o
rhs.o: params.o
