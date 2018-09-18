CC = g++ 
CR = gcc

CFLAGS = -g -std=c++11
CRFLAGS = 

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

.c.o:
	$(CR) $(CRFLAGS) -c $*.c -o $*.o

OBJS = gene.o model.o Rloop_equilibrium_model.o simulation.o structure.o windower.o

all: programs

programs: rlooper

bin:
	mkdir -p bin

rlooper: bin $(OBJS) ensemble_analyzer_main.o
	@echo Build rlooper.
	$(CC) -o bin/rlooper $(CFLAGS) $(OBJS) ensemble_analyzer_main.o
clean: 
	@echo Clean.
	rm -f *.o *.a tests/*.o

clobber: clean
	@echo Clobber.
	rm -fr bin

