CC = g++ 
CR = gcc

CFLAGS = -g -std=c++11
CRFLAGS = 

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

.c.o:
	$(CR) $(CRFLAGS) -c $*.c -o $*.o

OBJS = gene.o lumberjack.o model.o Rloop_equilibrium_model.o simulation.o structure.o windower.o

all: programs

programs: ensemble_analyzer perloop expected_length

bin:
	mkdir -p bin

ensemble_analyzer: bin $(OBJS) ensemble_analyzer_main.o
	@echo Build ensemble_analyzer.
	$(CC) -o bin/ensemble_analyzer $(CFLAGS) $(OBJS) ensemble_analyzer_main.o

perloop: bin $(OBJS) perloop_main.o
	@echo Build perloop.
	$(CC) -o bin/perloop  $(CFLAGS) $(OBJS) perloop_main.o

expected_length: bin $(OBJS) expected_length_main.o
	@echo Build expected_length
	$(CC) -o bin/expected_length $(CFLAGS) $(OBJS) expected_length_main.o

clean: 
	@echo Clean.
	rm -f *.o *.a tests/*.o xinger/*.o

clobber: clean
	@echo Clobber.
	rm -fr bin

