Timing your codes:

python examples:
       % ipython
       run -i timing1.py
       %time a=make1(1000,0)
       %time a=make1(1000,1)
       %time a=make1(1000,2)
       %time amake1(1000,3)
       %timeit make1(1000,3)

       %time d=delta1(a)
       %time d=delta2(a)

       %timeit d=delta1(a)
       %timeit d=delta2(a)


C example:
       code examples taken from NEMO package     http://carma.astro.umd.edu/nemo/

       # interesting (intel only) clock cycle routine in
       
       more $NEMO/src/kernel/misc/timers.c
       man timers




       cd $NEMO/src/nbody/evolve/hackcode/hackcode1
       make hackcode1 CC="gcc -pg"
       hackcode1 nbody=1000
       gprof hackcode1

       # -> not so great, but this is the classic  profiler

       # on mac, Xcode has a good profiler
       
       # intel compiler has a good profiler
       
       # linux has valgrid and kcachegrind
       valgrind --tool=callgrind  hackcode1 nbody=1000
       kcachegrind callgrind.out.16241

       # note valgrind has lots of other useful option (finding memory leaks etc.)
       # and also note the profiler doesn't need the -pg compile flag
       make hackcode1_qp
       valgrind --tool=callgrind  hackcode1_qp nbody=1000	
       kcachegrind callgrind.out.16301

       # note that the gravity routine now has twice as much CPU overhead.(25 -> 50%)

