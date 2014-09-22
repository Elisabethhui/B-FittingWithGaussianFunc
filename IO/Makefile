OBJS = IO.o OutputController.o
LIB=libIO.a

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
#FLAGS=-Wall -g  -I${HOME}/program/local/boost/include -L${HOME}/program/local/boost/lib
FLAGS=-Wall ${OPTCFLAGS}  -I${CPLUS_INCLUDE_PATH}
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
#	$(AR) -x $(SUBLIBS)
	$(AR) $(ARFLAGS) $@ *.o
		
cleanall:
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir cleanall || exit 1;\
	done
	$(RM) $(OBJS) $(LIB) *.bak *~

.PHONY: all cleanall
