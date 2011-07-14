LINESEARCH_OBJ = linesearch.o multilabel_linesearch.o hinge_linesearch.o

linesearch.o: ${LINESEARCH_DIR}/linesearch.hpp ${LINESEARCH_DIR}/linesearch.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${LINESEARCH_DIR}/linesearch.cpp

multilabel_linesearch.o: ${LINESEARCH_DIR}/multilabel_linesearch.hpp ${LINESEARCH_DIR}/multilabel_linesearch.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${LINESEARCH_DIR}/multilabel_linesearch.cpp

hinge_linesearch.o: ${LINESEARCH_DIR}/hinge_linesearch.hpp ${LINESEARCH_DIR}/hinge_linesearch.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${LINESEARCH_DIR}/hinge_linesearch.cpp
