
  COMP   = $(MPI)/bin/mpic++
  EXT    = cpp

SRC_PATH= ../src/*.cpp
OUT_PATH= ../bin/
OBJ_PATH= ../obj/
INC_PATH= ../include/
EXE_PATH= ../../bin/

  MPI  = /usr/local

 PROG = GITR
 
##--------------------> set code components: 
  CODE_o = z.main.o z.mainMR.o z.mainWR.o z.io.o z.setup.o z.messaging.o Particle.o EfieldInterp.o

##-------------------- not needed for most 'make':
# CODE  = z.main.$(EXT) z.mainMR.$(EXT) z.mainWR.$(EXT) z.io.$(EXT)  \
#         z.setup.$(EXT) z.update.$(EXT) z.messaging.$(EXT)
 
#======================= create executable: make compile ============# 


compile: 
	$(COMP) $(FLAGs) -I$(INC_PATH) $(SRC_PATH)  -o $(OUT_PATH)$(PROG).x 

#======================= execute: make run | make pbs ==========# 
run: 
	$(MPI)/bin/mpirun -n  5  $(EXE_PATH)$(PROG).x

