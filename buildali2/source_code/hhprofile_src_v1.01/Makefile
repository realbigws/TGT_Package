CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -fno-strict-aliasing -static
OBJ = hhpred_ali.o hhpred_hit.o hhpred_hmm.o hhpred_mat.o hhpred_util.o

all: hh_filter hh_make
	rm -f $(OBJ) hh_filter.o hh_make.o

hh_filter: $(OBJ) hh_filter.o
	${CXX} ${CCFLAGS} -o hh_filter $(OBJ) hh_filter.o

hh_make: $(OBJ) hh_make.o
	${CXX} ${CCFLAGS} -c hh_make.cpp
	${CXX} ${CCFLAGS} -o hh_make $(OBJ) hh_make.o


hh_filter.o : 
	${CXX} ${CCFLAGS} -c hh_filter.cpp
hh_make.o :
	${CXX} ${CCFLAGS} -c hh_make.cpp
hhpred_ali.o :
	${CXX} ${CCFLAGS} -c hhpred_ali.cpp
hhpred_hit.o :
	${CXX} ${CCFLAGS} -c hhpred_hit.cpp
hhpred_hmm.o :
	${CXX} ${CCFLAGS} -c hhpred_hmm.cpp
hhpred_mat.o :
	${CXX} ${CCFLAGS} -c hhpred_mat.cpp
hhpred_util.o :
	${CXX} ${CCFLAGS} -c hhpred_util.cpp
clean:
	rm -f hh_filter hh_make $(OBJ) hh_filter.o hh_make.o

