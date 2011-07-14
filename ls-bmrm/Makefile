BOOST_DIR=../externalpackages/
UTILITIES_DIR=../utilities
DATA_DIR=../data
MODEL_DIR=../model
LOSS_DIR=../loss
SOLVER_DIR=../solver
INNER_SOLVER_DIR=../solver/bmrminnersolver
LINESEARCH_DIR=../linesearch

include ${LOSS_DIR}/linesearchloss.mk
include ${DATA_DIR}/data.mk
include ${LINESEARCH_DIR}/linesearch.mk
include ${SOLVER_DIR}/solver.mk
include ${UTILITIES_DIR}/utilities.mk
include ${MODEL_DIR}/model.mk

FACTORY_SRC = ${LOSS_DIR}/linesearchlossfactory.hpp ${DATA_DIR}/datafactory.hpp\
	      ${SOLVER_DIR}/solverfactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp


INCLUDE_DIR = -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR}\
	      -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR} -I${LINESEARCH_DIR}

ALL_OBJECTS = ${UTILITIES_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LINESEARCH_LOSS_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${LINESEARCH_OBJ}

CXX = g++

### for debugging ###
#CFLAGS = -Wall -O0 -g -DDEBUG ${INCLUDE_DIR}

### for release ###
CFLAGS = -Wall -O3 -DNDEBUG ${INCLUDE_DIR} -DLINESEARCH_BMRM

LDFLAGS = 


.PHONY: all

all: ls-bmrm-train ls-bmrm-predict

# executable
ls-bmrm-train: ls-bmrm-train.cpp ${ALL_OBJECTS}  ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} ls-bmrm-train.cpp -o ls-bmrm-train ${ALL_OBJECTS}

ls-bmrm-predict: ls-bmrm-predict.cpp ${ALL_OBJECTS} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} ls-bmrm-predict.cpp -o ls-bmrm-predict ${ALL_OBJECTS}


clean:
	rm *.o ls-bmrm-train ls-bmrm-predict;


