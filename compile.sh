ifort -Vaxlib  -extend_source -g -heap-arrays -traceback -mcmodel=large -i-dynamic -g -r8 -i8 -O2    -c galtogad.f90 makehalo.f90 routines.f90
ifort -Vaxlib  -extend_source -traceback -fpe0  -mcmodel=large -i-dynamic -g -r8 -i8 -O2 galtogad.o  makehalo.o routines.o -o galtogadrun
chmod u+x galtogadrun
#ifort -Vaxlib  -extend_source -g -heap-arrays -traceback  -mcmodel=large -i-dynamic -r8 -i8 -O2 -c galtogad_simpletemp.f90  makehalo.f90  routines.f90
#ifort -Vaxlib  -extend_source  -traceback -fpe0  -mcmodel=large -i-dynamic -r8 -i8 -O2 galtogad_simpletemp.o  makehalo.o  routines.o -o galtogadrun_simpletemp
#chmod u+x galtogadrun
