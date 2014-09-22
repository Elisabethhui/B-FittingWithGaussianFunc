TARGET:= Fitting.exe
OBJS := main.o
LIB:=IO.a

SUBDIRS:=Math  GaussFitting Common Constants IO Parser Parallelization

SUBLIBS:=\
		 ./GaussFitting/libGaussFitting.a \
		 ./Math/libMath.a \
		./IO/libIO.a \
		./Parser/libParser.a \
		./Parallelization/libParallelization.a \
		./Common/libCommon.a \
		./Constants/libConstants.a \

RM := -rm -rf
CC=mpicc

ifeq (${DEBUG}, y)
OPTCFLAGS += -g
else
OPTCFLAGS += -O2
endif

ifeq (${PROFILEFLAGS},y)
OPTCFLAGS += -pg
endif


CFLAGS=-Wall ${OPTCFLAGS}  -I${HOME}/program/local/gsl/include

AR:=ar
ARFLAGS:=crvs

all:${TARGET}

$(TARGET):${OBJS}
	for dir in $(SUBDIRS);\
	do \
		make  -C $$dir all || exit 1;\
	done
	$(CC) -o main.o -c ${CFLAGS} main.cpp
	${CC}  -o $@ main.o ${CFLAGS} ${SUBLIBS} -L${HOME}/program/local/boost/lib -lboost_system -lboost_filesystem -L/${HOME}/program/local/gsl/lib -lgsl -lfftw3 -lgslcblas -lm

		
cleanall:
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir cleanall || exit 1;\
	done
	$(RM) $(OBJS) $(LIB) *.bak *~ *.o ${TARGET}

.PHONY: all cleanall
