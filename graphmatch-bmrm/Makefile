# Set to 1 if COIN Clp is installed and you would like to use L1 norm regularizer
COIN_CLP_AVAILABLE=0

# Set to 1 if BT code by Schramm and Zowe (1992) is avaialble and you would like to use it
BT_AVAILABLE=0


BOOST_DIR=/usr/include
UTILITIES_DIR=../utilities
DATA_DIR=../data
MODEL_DIR=../model
LOSS_DIR=../loss
SOLVER_DIR=../solver
INNER_SOLVER_DIR=../solver/bmrminnersolver
ASP_DIR=asp


include ${SOLVER_DIR}/solver.mk
include ${UTILITIES_DIR}/utilities.mk
include ${MODEL_DIR}/model.mk

DATA_OBJ = data.o graphdata.o
LOSS_OBJ = loss.o graphmatchloss.o
ASP_OBJ = lap.o gap.o

FACTORY_SRC = ${SOLVER_DIR}/solverfactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp


INCLUDE_DIR = -I${ASP_DIR} -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR}

CXX = g++

ALL_OBJECTS = ${UTILITIES_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${ASP_OBJ}

### for debugging ###
#CFLAGS = -Wall -O0 -g -DNDEBUG  ${INCLUDE_DIR}

### for release ###
CFLAGS = -Wall -O3 -DNDEBUG ${INCLUDE_DIR}

ifeq (${COIN_CLP_AVAILABLE},1)
	CFLAGS += ${COIN_CFLAGS} -DHAVE_L1N1_INNER_SOLVER
	LDFLAGS += ${COIN_LDFLAGS}
endif

ifeq (${BT_AVAILABLE},1)
	REGULARIZER_DIR=../regularizer
	include ${REGULARIZER_DIR}/regularizer.mk
	FACTORY_SRC += ${REGULARIZER_DIR}/regularizerfactory.hpp
	INCLUDE_DIR += -I${REGULARIZER_DIR}
	ALL_OBJECTS += ${REGULARIZER_OBJ}
	CFLAGS += -DHAVE_BT_SOLVER
	LDFLAGS += ${BT_LDFLAGS}
endif

.PHONY: all

all: gm-bmrm-train gm-bmrm-predict

# executable
gm-bmrm-train: gm-bmrm-train.cpp ${ALL_OBJECTS}  ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} gm-bmrm-train.cpp -o gm-bmrm-train ${ALL_OBJECTS} CImg.h -lm -lpthread -lX11

gm-bmrm-predict: gm-bmrm-predict.cpp ${ALL_OBJECTS} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} gm-bmrm-predict.cpp -o gm-bmrm-predict ${ALL_OBJECTS} CImg.h -lm -lpthread -lX11



# loss objects
loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

graphmatchloss.o: graphmatchloss.hpp graphmatchloss.cpp graphdata.cpp gap.o lap.o
	${CXX} ${CFLAGS} -c graphmatchloss.cpp


# Assignment solvers.
lap.o: ${ASP_DIR}/lap.cpp ${ASP_DIR}/lap.hpp
	${CXX} ${CFLAGS} -c ${ASP_DIR}/lap.cpp

gap.o: ${ASP_DIR}/gap.h ${ASP_DIR}/gap.c
	${CXX} ${CFLAGS} -c ${ASP_DIR}/gap.c


# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

graphdata.o: graphdata.hpp graphdata.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c graphdata.cpp


clean:
	rm *.o gm-bmrm-train gm-bmrm-predict;


