REGULARIZER_OBJ = regularizer.o l1n1.o l2n2.o

# regularizer objects
regularizer.o: ${REGULARIZER_DIR}/regularizer.hpp ${REGULARIZER_DIR}/regularizer.cpp
	${CXX} ${CFLAGS} -c ${REGULARIZER_DIR}/regularizer.cpp

l1n1.o: ${REGULARIZER_DIR}/l1n1.hpp ${REGULARIZER_DIR}/l1n1.cpp
	${CXX} ${CFLAGS} -c ${REGULARIZER_DIR}/l1n1.cpp

l2n2.o: ${REGULARIZER_DIR}/l2n2.hpp ${REGULARIZER_DIR}/l2n2.cpp
	${CXX} ${CFLAGS} -c ${REGULARIZER_DIR}/l2n2.cpp
