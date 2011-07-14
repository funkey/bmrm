UTILITIES_OBJ = sml.o common.o configuration.o bmrmexception.o timer.o

# utility objects
sml.o: ${UTILITIES_DIR}/sml.hpp ${UTILITIES_DIR}/sml.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/sml.cpp

common.o: ${UTILITIES_DIR}/common.hpp ${UTILITIES_DIR}/common.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/common.cpp

configuration.o: ${UTILITIES_DIR}/configuration.hpp ${UTILITIES_DIR}/configuration.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/configuration.cpp

bmrmexception.o: ${UTILITIES_DIR}/bmrmexception.hpp ${UTILITIES_DIR}/bmrmexception.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/bmrmexception.cpp

timer.o: ${UTILITIES_DIR}/timer.hpp ${UTILITIES_DIR}/timer.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/timer.cpp
