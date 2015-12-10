.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -fpp $(FOPT) $(INC) $(IOTK_INC) -c $<

F90=ifort

FOPT=-assume byterecl

ESPRESSO_DIR=/public/home/ccao/source/espresso-5.2.1
IOTK_DIR=$(ESPRESSO_DIR)/iotk

IOTK_INC=-I$(IOTK_DIR)/src
IOTK_LIB=$(IOTK_DIR)/src/libiotk.a

#LIB=-L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

all: qe vasp

qe : z2inv.f90 wave_pw.o constants.o
	$(F90) $(FOPT) -o z2inv_qe.x z2inv.f90 wave_pw.o constants.o $(IOTK_INC) $(IOTK_LIB)

vasp : z2inv.f90 wave_vasp.o constants.o
	$(F90) $(FOPT) -o z2inv_vasp.x z2inv.f90 wave_vasp.o constants.o

wave_pw.o: wave_pw.f90 constants.o

wave_vasp.o: wave_vasp.f90 constants.o

constants.o: constants.f90

clean:
	rm -rf *.x *.mod *.o
