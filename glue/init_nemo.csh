#! /bin/csh -f

#   generate data using NEMO (http://www.astro.umd.edu/nemo)

set nb=100000
set ns=128

echo x y z vx vy vz > model1.tab
mkdisk - 1000 z0=0.5 seed=123 rmax=2 | snapprint - >> model1.tab

echo x y vx vy > model2.tab
mkspiral - mass=1 a=100 w=5 seed=123 rmax=2 |\
 snapprint - x,y,vx,vy >> model2.tab

rm -f model2.fits
mkspiral - $nb mass=1 a=100 w=5 seed=123 rmax=2 |\
 snapgrid - - nx=$ns ny=$ns                              |\
  ccdfits - model2i.fits radecvel=true

rm -f model3.fits
mkspiral - $nb mass=1 a=100 w=5 seed=123 rmax=2 |\
 snapgrid - - nx=$ns ny=$ns nz=$ns/2 zrange=-1:1 zvar=vy |\
  ccdfits - model3i.fits radecvel=true
