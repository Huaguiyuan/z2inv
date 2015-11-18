F90=ifort

ESPRESSO_DIR=/public/home/ccao/source/espresso-5.2.1
IOTK_DIR=$(ESPRESSO_DIR)/iotk

INC=-I$(IOTK_DIR)/src -I$(ESPRESSO_DIR)/Modules
LIB=$(IOTK_DIR)/src/libiotk.a

z2inv.x : z2inv.f90
	$(F90) -o z2inv.x z2inv.f90 $(INC) $(LIB)
