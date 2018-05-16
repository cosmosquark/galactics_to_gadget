OPTIONS = -O -DRINGASCII -DASCII -m32
SYS = F
CPPFLAGS = $(OPTIONS) -D$(SYS)
#FLAGS =  -autodouble -openmp -Vaxlib -extend_source -mcmodel=large -i-dynamic  
FLAGS = -Vaxlib -openmp -extend_source -mcmodel=large -i-dynamic -g -r8 -i8
#FLAGS = -O2 -mcmodel=medium
#FLAGS = -O2 -mcmodel=medium -shared-intel
F77=ifort
#F77 = mpif77

.SUFFIXES: .F .F90

.F90.o:
	$(F77) $(FLAGS) -c $<
.f.o:
	$(F77) $(FLAGS) -c -openmp $<
.F.o:
	@ if test $(SYS) = ibm ; then \
	echo "$(CPP) -P $(CPPFLAGS) $< > ../tmp/$*.f";\
	$(CPP) -P $(CPPFLAGS) $< > ../tmp/$*.f;\
	echo "$(F77) $(FLAGS) -c ../tmp/$*.f";\
	$(F77) $(FLAGS) -c ../tmp/$*.f;\
        elif test $(SYS) = linux ; then \
		echo "$(CPP) -P $(CPPFLAGS) $< > ../tmp/$*.f";\
		$(CPP) -P $(CPPFLAGS) $< > ../tmp/$*.f;\
		echo "$(F77) $(FLAGS) -c ../tmp/$*.f";\
		$(F77) $(FLAGS) -c ../tmp/$*.f;\
	else \
	echo "$(F77) $(FLAGS) $(CPPFLAGS) -c $<";\
	$(F77) $(FLAGS) $(CPPFLAGS) -c $<;\
	fi

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

#
galacticstogadget: galacticstogadget.o ahtime.o gravsubs.o dumpdata2sf.o system.o inunit.o createcool.o sortrx.o kernel.o makehalo.o calculate_h0_fun.o ssort.o
	$(F77) $(FLAGS) twogalaxiesmhh.o tiltgal_mhh.o ahtime.o gravsubs.o dumpdata2sf.o system.o inunit.o createcool.o sortrx.o kernel.o makehalo.o calculate_h0_fun.o ssort.o -o twogalaxiesmhh
	mv twogalaxiesmhh ../$(RUNDIR)

clean:
	/bin/rm  -rf *.o tmp

system.o: system.$(SYS)
	@ if test $(SYS) = linux ; then \
	cp system.linux system.c; \
	$(CC) -c $(CPPFLAGS) $(CFLAGS) system.c; \
	else \
	cp system.$(SYS) system.f; \
	$(F77) $(FLAGS) -c system.f; \
	fi

