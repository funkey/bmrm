# Set to 1 if COIN Clp is installed
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

include ${UTILITIES_DIR}/utilities.mk
include ${MODEL_DIR}/model.mk
include ${SOLVER_DIR}/solver.mk

DATA_OBJ = data.o seqmulticlassfeature.o seqmulticlasslabel.o seqmulticlassdata.o
LOSS_OBJ = loss.o smmmulticlassloss.o

FACTORY_SRC = ${SOLVER_DIR}/solverfactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp

INCLUDE_DIR =  -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR}

ALL_OBJECTS = ${UTILITIES_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ}

CXX = g++


### for debugging ###
#CFLAGS = -O0 -g  -Wall -DNDEBUG

### for release ###
CFLAGS = -O3 -Wall -DNDEBUG ${INCLUDE_DIR}

LDFLAGS = 

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

all: smmmc-bmrm-train smmmc-bmrm-predict

# executable
smmmc-bmrm-train: smmmc-bmrm-train.cpp ${ALL_OBJECTS} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} smmmc-bmrm-train.cpp -o smmmc-bmrm-train ${ALL_OBJECTS}

smmmc-bmrm-predict: smmmc-bmrm-predict.cpp ${ALL_OBJECTS} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} smmmc-bmrm-predict.cpp -o smmmc-bmrm-predict ${ALL_OBJECTS} 


# loss objects
loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

smmmulticlassloss.o: smmmulticlassloss.hpp smmmulticlassloss.cpp seqmulticlassdata.o
	${CXX} ${CFLAGS} -c smmmulticlassloss.cpp


# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

seqmulticlassfeature.o: seqmulticlassfeature.hpp seqmulticlassfeature.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c seqmulticlassfeature.cpp

seqmulticlasslabel.o: seqmulticlasslabel.hpp seqmulticlasslabel.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c seqmulticlasslabel.cpp

seqmulticlassdata.o: seqmulticlassdata.hpp seqmulticlassdata.cpp  seqmulticlassfeature.o seqmulticlasslabel.o ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c seqmulticlassdata.cpp


clean:
	rm *.o smmmc-bmrm-train smmmc-bmrm-predict;


