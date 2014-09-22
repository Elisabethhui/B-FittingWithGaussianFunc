OBJS := Header.o
LIB :=libCommon.a

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
FLAGS=-Wall ${OPTCFLAGS}  -lfftw3 -lm
AR=ar
ARFLAGS=crvs


all: $(LIB)

$(OBJS):%.o:%.cpp
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir all || exit 1;\
	done
	$(CC) -c $(FLAGS) -o $@  $<

$(LIB):$(OBJS)
	$(AR) $(ARFLAGS) $@ *.o
		
cleanall:
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir cleanall || exit 1;\
	done
	$(RM) $(OBJS) $(LIB) *.bak *~ *.o

.PHONY: all cleanall
