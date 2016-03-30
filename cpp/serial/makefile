
COMP   = g++

SRC_PATH= src
OUT_PATH= bin
OBJ_PATH= obj
INC_PATH= include
#path for libconfig
LIB_PATH= ~/libconfig-1.5/lib/

PROG = GITR_serial
 
compile: 
	$(COMP) $(FLAGs) -I$(INC_PATH) -I$(LIB_PATH) $(SRC_PATH)/*.cpp  -o $(OUT_PATH)/$(PROG).x -lconfig++

clean: 
	rm $(OUT_PATH)/* $(OBJ_PATH)/* $(OUT_PATH)/*	

