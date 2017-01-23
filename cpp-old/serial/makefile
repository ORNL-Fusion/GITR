#for local compile on mac
#g++  -Iinclude -I~/libconfig-1.5/lib/ src/*.cpp  -o bin/GITR_serial.x -lconfig++ -std=c++11
COMP   = g++

SRC_PATH= src
OUT_PATH= bin
OBJ_PATH= obj
INC_PATH= include
#path for libconfig
#LIB_PATH= ~/libconfig-1.5/lib/
LIB_PATH= /home/dg6/code/libconfig/gnu/lib
PROG = GITR_serial
#INC_LIB = /home/dg6/code/libconfig/gnu/include
 
compile: 
	$(COMP) $(FLAGs) -I$(INC_PATH) -L$(LIB_PATH) $(SRC_PATH)/*.cpp  -o $(OUT_PATH)/$(PROG).x -I$(INC_LIB) -lconfig++ -std=c++11

clean: 
	rm $(OUT_PATH)/* $(OBJ_PATH)/* $(OUT_PATH)/*	

