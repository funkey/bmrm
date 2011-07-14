DATA_OBJ = data.o vecfeature.o veclabel.o varlenveclabel.o vecdata.o multilabelvecdata.o

# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

vecfeature.o: ${DATA_DIR}/vecfeature.hpp ${DATA_DIR}/vecfeature.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecfeature.cpp

veclabel.o: ${DATA_DIR}/veclabel.hpp ${DATA_DIR}/veclabel.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/veclabel.cpp

varlenveclabel.o: ${DATA_DIR}/varlenveclabel.hpp ${DATA_DIR}/varlenveclabel.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/varlenveclabel.cpp

vecdata.o: ${DATA_DIR}/vecdata.hpp ${DATA_DIR}/vecdata.cpp data.o vecfeature.o veclabel.o ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecdata.cpp

multilabelvecdata.o: ${DATA_DIR}/multilabelvecdata.hpp ${DATA_DIR}/multilabelvecdata.cpp data.o vecfeature.o varlenveclabel.o ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/multilabelvecdata.cpp
