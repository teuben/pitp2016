#

FC  = g77
CC  = gcc
CXX = g++
FLIBS = -lg2c


ARCH = `uname -s`_`uname -m`

-include makedefs

BINS = hello1 hello2 hello3 hello12_c hello12_f hello21_c hello21_f

all: $(BINS)

subs:


hello1:
	$(CC) hello1.c    -o hello1

hello2:
	$(FC) hello2.f    -o hello2

hello3:
	$(CXX) hello3.cc  -o hello3


hello12_c:
	$(CC) -c main1.c
	$(FC) -c sub2.f
	$(CC) main1.o sub2.o   $(FLIBS)  -o hello12_c

hello12_f:
	$(CC) -c main1.c
	$(FC) -c sub2.f
	$(FC) main1.o sub2.o -o hello12_f

hello21_c:
	$(FC) -c main2.f
	$(CC) -c sub1.c
	$(CC) main2.o sub1.o  $(FLIBS) -o hello21_c

hello21_f:
	$(FC) -c main2.f
	$(CC) -c sub1.c
	$(FC) main2.o sub1.o -o hello21_f




tests: $(BINS)
	-@for b in $(BINS); do \
	 ./$$b ; done

clean:
	rm -rf $(BINS) *.o


save:
	mkdir -p $(ARCH)
	cp $(BINS) $(ARCH)
	uname -a > $(ARCH)/uname.log

test2:
	-@for b in $(BINS); do \
	 $(ARCH)/$$b ; done
	
