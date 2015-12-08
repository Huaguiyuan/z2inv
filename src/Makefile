.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -fpp $(FOPT) $(INC) -c $<

F90=ifort

FOPT=-assume byterecl
INC=-I.

#ESPRESSO_DIR=/public/home/ccao/source/espresso-5.2.1
#IOTK_DIR=$(ESPRESSO_DIR)/iotk

#INC=-I$(IOTK_DIR)/src -I$(ESPRESSO_DIR)/Modules
#LIB=$(IOTK_DIR)/src/libiotk.a

#z2inv.x : z2inv.f90
#	$(F90) -o z2inv.x z2inv.f90 $(INC) $(LIB)

z2inv_vasp.x : z2inv_vasp.f90 wvfunc.o constants.o
	$(F90) -fpp $(FOPT) $(INC) -o z2inv_vasp.x z2inv_vasp.f90 wvfunc.o constants.o flib.o

constants.o: constants.f90
flib.o: flib.f90

wvfunc.o: wvfunc.f90 constants.o flib.o


clean:
	rm -rf *.x *.mod *.o
