OBJS =Math.o sphere_lebedev_rule.o
LIB=libMath.a

ifeq ($(DEBUG), y)
OPTCFLAGS += -g
else
OPTCFLAGS += -O2
endif

ifeq (${PROFILEFLAGS},y)
OPTCFLAGS += -pg
endif



SUBDIRS:=
SUBLIBS:=
RM = -rm -rf
CC=mpicc
CFLAGS=-Wall ${OPTCFLAGS} -I/usr/local/fftw/include -L/usr/local/fftw/lib -lfftw3 -I${HOME}/program/local/gsl/include -L${HOME}/program/local/gsl/lib -lgsl -lgslcblas -lm 
AR=ar
ARFLAGS=crvs


all: $(OBJS) $(LIB)

$(OBJS):%.o:%.cpp
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir all || exit 1;\
	done
	${CC} -c ${CFLAGS}  -o $@  $< 

$(LIB):$(OBJS)
#	$(AR) -x $(SUBLIBS)
	$(AR) $(ARFLAGS) $@ *.o
		
cleanall:
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir cleanall || exit 1;\
	done
	$(RM) $(OBJS) $(LIB) *.bak *~ *.o

.PHONY: all cleanall
