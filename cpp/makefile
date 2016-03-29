
COMP   = nvcc

SRC_PATH= src
OUT_PATH= bin
OBJ_PATH= obj
INC_PATH= include
#path for libconfig header
HED_PATH= /usr/include/
#path for libconfig++.so
LIB_PATH= /home/dg6/code/libconfig/source/libconfig-1.5/lib/.libs/
#thrust header files
#THR_PATH = /cm/shared/apps/cuda70/toolkit/7.0.28/include
THR_PATH = /cm/shared/apps/cuda75/toolkit/7.5.18/include
PROG = GITR_serial
 
compile: 
	$(COMP) $(FLAGs) -I$(INC_PATH) -I$(HED_PATH) -I$(THR_PATH) -L$(LIB_PATH) -lconfig++ $(SRC_PATH)/*.cpp  -o $(OUT_PATH)/$(PROG).x

clean: 
	rm $(OUT_PATH)/* $(OBJ_PATH)/* $(OUT_PATH)/*	

