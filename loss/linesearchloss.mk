LINESEARCH_LOSS_OBJ = loss.o scalarloss.o binaryclassificationloss.o linesearch_hingeloss.o\
						multilabelloss.o linesearch_multilabelloss.o

loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

scalarloss.o: ${LOSS_DIR}/scalarloss.hpp ${LOSS_DIR}/scalarloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/scalarloss.cpp

binaryclassificationloss.o: ${LOSS_DIR}/binaryclassificationloss.hpp ${LOSS_DIR}/binaryclassificationloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/binaryclassificationloss.cpp

multilabelloss.o: ${LOSS_DIR}/multilabelloss.cpp ${LOSS_DIR}/multilabelloss.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/multilabelloss.cpp

linesearch_hingeloss.o: ${LOSS_DIR}/linesearch_hingeloss.hpp ${LOSS_DIR}/linesearch_hingeloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/linesearch_hingeloss.cpp

linesearch_multilabelloss.o: ${LOSS_DIR}/linesearch_multilabelloss.cpp ${LOSS_DIR}/linesearch_multilabelloss.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/linesearch_multilabelloss.cpp
