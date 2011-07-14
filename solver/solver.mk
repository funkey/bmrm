SOLVER_OBJ = bmrm.o

# NOTE:
# gfortran is required to compile the original bt fortran code
#
ifeq (${BT_AVAILABLE},1)
bt.o: ${SOLVER_DIR}/bt.hpp ${SOLVER_DIR}/bt.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/bt.cpp  ${DEBUGFLAG}

bt_orig.o: ${SOLVER_DIR}/bt_orig.f
	gfortran -c -O3 ${SOLVER_DIR}/bt_orig.f

BT_LDFLAGS = -lgfortran
SOLVER_OBJ += bt.o bt_orig.o
endif


INNER_SOLVER_OBJ = bmrminnersolver.o l2n2_bmrmdualinnersolver.o\
	           l2n2_daifletcherpgm.o l2n2_prloqo.o l2n2_linesearch.o

# NOTE: 
# Environment variables COIN_INC_DIR and COIN_LIB_DIR must be set to the
# include/ and lib/ folders under the COIN Clp installation, respectively.
#
ifeq (${COIN_CLP_AVAILABLE},1)
COIN_CFLAGS = -I${COIN_INC_DIR} -O3 -fomit-frame-pointer -pipe -DNDEBUG  -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion 

COIN_LDFLAGS = -L${COIN_LIB_DIR} -lm -lClp -lCoinUtils `cat ${COIN_LIB_DIR}/coinutils_addlibs.txt`

l1n1_clp.o: ${INNER_SOLVER_DIR}/l1n1_clp.hpp ${INNER_SOLVER_DIR}/l1n1_clp.cpp
	${CXX} ${CFLAGS} ${COIN_CFLAGS} -c ${INNER_SOLVER_DIR}/l1n1_clp.cpp

INNER_SOLVER_OBJ += l1n1_clp.o
endif



# solver objects
bmrm.o: ${SOLVER_DIR}/bmrm.hpp ${SOLVER_DIR}/bmrm.cpp
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/bmrm.cpp

bmrminnersolver.o: ${INNER_SOLVER_DIR}/bmrminnersolver.hpp ${INNER_SOLVER_DIR}/bmrminnersolver.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/bmrminnersolver.cpp

l2n2_bmrmdualinnersolver.o: ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.hpp ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp bmrminnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp

l2n2_daifletcherpgm.o: ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.hpp ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp

l2n2_prloqo.o: ${INNER_SOLVER_DIR}/l2n2_prloqo.hpp ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp 

l2n2_linesearch.o: ${INNER_SOLVER_DIR}/l2n2_linesearch.hpp ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp 


