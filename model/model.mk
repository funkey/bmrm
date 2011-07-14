MODEL_OBJ = model.o

# model objects
model.o: ${MODEL_DIR}/model.hpp ${MODEL_DIR}/model.cpp ${UTILITIES_OBJ}
	${CXX} ${CFLAGS} -c ${MODEL_DIR}/model.cpp

