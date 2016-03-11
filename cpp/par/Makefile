
#============================= set MPI, compiler ======================#
##-----> set appropriate compiler_wrapper: mpif77 mpif90 mpicc mpic++
  COMP   = $(MPI)/bin/mpic++
##-----> set appropriate extension: f  c  cpp
  EXT    = cpp
  LFLAGs = 
#for C:  LFLAGs = -lm
##-------------------------- for all:
  FLAGs  = -g  $(MPIlnk) 
# FLAGs  = -O3 $(MPIlnk) 
  MPIlnk = -I$(MPI)/include  -L$(MPI)/lib
##---------------------> set path to openMPI on local:
  MPI  = /usr/local
##---------------------> parallel on cluster:  for open-mpi with intel
# MPI  = /opt/open-mpi/tcp-intel11
##-------------------- parallel on cluster:  for open-mpi with gnu44
# MPI  = /opt/open-mpi/tcp-gnu44
#
#=========================== set source code  =========================#
##--------------->set names for your PROGram and std I/O files: 
#  PROG = code1D
 PROG = code2D
  INPUT  = ./dat.txt 
  OUTPUT = o.out 
##--------------------> set code components: 
  CODE_o = z.main.o z.mainMR.o z.mainWR.o z.io.o z.setup.o z.messaging.o Particle.o EfieldInterp.o

##-------------------- not needed for most 'make':
# CODE  = z.main.$(EXT) z.mainMR.$(EXT) z.mainWR.$(EXT) z.io.$(EXT)  \
#         z.setup.$(EXT) z.update.$(EXT) z.messaging.$(EXT)
 
#======================= create executable: make compile ============# 
$(CODE_o):%.o: %.$(EXT)
	$(COMP) $(FLAGs) -c $< -o $@

compile: $(CODE_o) 
	$(COMP) $(FLAGs)  $(CODE_o)  -o $(PROG).x  $(LFLAGs)
	@echo " >>> compiled on `hostname -s` with  $(COMP) <<<" 
	rm -f z.*.o 
	rm Particle.o
#======================= execute: make run | make pbs ==========# 
run: 
	$(MPI)/bin/mpirun -n  5  $(PROG).x
#	 < $(INPUT) > $(OUTPUT) 
#or:	$(MPI)/bin/mpirun  -np 3  $(PROG).x < $(INPUT) > $(OUTPUT) 
pbs: 
	@ vi PBSscript 
	make clean 
	qsub  PBSscript 

clean: 
	rm -f  o.*  DONE  watch 
#----------------------------------------------------------------------#
