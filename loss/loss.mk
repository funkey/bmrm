LOSS_OBJ = loss.o scalarloss.o\
	   binaryclassificationloss.o hingeloss.o squaredhingeloss.o logisticloss.o rocscoreloss.o fbetaloss.o exponentialloss.o huberhingeloss.o\
	   univariateregressionloss.o leastsquaresloss.o epsiloninsensitiveloss.o quantileloss.o leastabsdevloss.o huberrobustloss.o\
	   noveltyloss.o poissonloss.o\
	   wtamulticlassloss.o multilabelloss.o\
	   rankloss.o ndcgrankloss.o lap.o

# loss objects
loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

scalarloss.o: ${LOSS_DIR}/scalarloss.hpp ${LOSS_DIR}/scalarloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/scalarloss.cpp

binaryclassificationloss.o: ${LOSS_DIR}/binaryclassificationloss.hpp ${LOSS_DIR}/binaryclassificationloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/binaryclassificationloss.cpp

univariateregressionloss.o: ${LOSS_DIR}/univariateregressionloss.hpp ${LOSS_DIR}/univariateregressionloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/univariateregressionloss.cpp

rankloss.o: ${LOSS_DIR}/rankloss.hpp ${LOSS_DIR}/rankloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rankloss.cpp

hingeloss.o: ${LOSS_DIR}/hingeloss.hpp ${LOSS_DIR}/hingeloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/hingeloss.cpp

squaredhingeloss.o: ${LOSS_DIR}/squaredhingeloss.hpp ${LOSS_DIR}/squaredhingeloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/squaredhingeloss.cpp

huberhingeloss.o: ${LOSS_DIR}/huberhingeloss.hpp ${LOSS_DIR}/huberhingeloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/huberhingeloss.cpp

logisticloss.o: ${LOSS_DIR}/logisticloss.hpp ${LOSS_DIR}/logisticloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/logisticloss.cpp

rocscoreloss.o: ${LOSS_DIR}/rocscoreloss.hpp ${LOSS_DIR}/rocscoreloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rocscoreloss.cpp

fbetaloss.o: ${LOSS_DIR}/fbetaloss.hpp ${LOSS_DIR}/fbetaloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/fbetaloss.cpp

exponentialloss.o: ${LOSS_DIR}/exponentialloss.hpp ${LOSS_DIR}/exponentialloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/exponentialloss.cpp

epsiloninsensitiveloss.o: ${LOSS_DIR}/epsiloninsensitiveloss.hpp ${LOSS_DIR}/epsiloninsensitiveloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/epsiloninsensitiveloss.cpp

leastsquaresloss.o: ${LOSS_DIR}/leastsquaresloss.hpp ${LOSS_DIR}/leastsquaresloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/leastsquaresloss.cpp

leastabsdevloss.o: ${LOSS_DIR}/leastabsdevloss.hpp ${LOSS_DIR}/leastabsdevloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/leastabsdevloss.cpp

huberrobustloss.o: ${LOSS_DIR}/huberrobustloss.hpp ${LOSS_DIR}/huberrobustloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/huberrobustloss.cpp

quantileloss.o: ${LOSS_DIR}/quantileloss.hpp ${LOSS_DIR}/quantileloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/quantileloss.cpp

noveltyloss.o: ${LOSS_DIR}/noveltyloss.hpp ${LOSS_DIR}/noveltyloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/noveltyloss.cpp

poissonloss.o: ${LOSS_DIR}/poissonloss.hpp ${LOSS_DIR}/poissonloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/poissonloss.cpp

wtamulticlassloss.o: ${LOSS_DIR}/wtamulticlassloss.hpp ${LOSS_DIR}/wtamulticlassloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/wtamulticlassloss.cpp

multilabelloss.o: ${LOSS_DIR}/multilabelloss.cpp ${LOSS_DIR}/multilabelloss.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/multilabelloss.cpp

ndcgrankloss.o: ${LOSS_DIR}/ndcgrankloss.hpp ${LOSS_DIR}/ndcgrankloss.cpp lap.o 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/ndcgrankloss.cpp

lap.o: ${LOSS_DIR}/lap.cpp ${LOSS_DIR}/lap.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/lap.cpp
