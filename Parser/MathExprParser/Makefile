OBJS =  error.o functions.o parser.o variablelist.o	 
RM = -rm -rf
CC=mpicc
FLAGS=-Wall -g
LIB=parser.a
AR=ar
ARFLAGS=cr

all: $(LIB)

$(LIB):$(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
		
cleanall:
	for dir in $(SUBDIRS);\
	do \
		make -C $$dir cleanall || exit 1;\
	done
	$(RM) $(OBJS) $(LIB) *.bak *~

.PHONY: all cleanall
