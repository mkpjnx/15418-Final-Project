CXX=g++ -m64
CXXFLAGS=-Icommon -Iobjs/ -O3 -Wall -g
MPI=-DMPI
MPICC = mpicxx

APP_NAME=grayscott
OBJDIR=objs

CXXFILES = main.cpp sim.cpp setup.cpp instrument.cpp
DEPS = cycleTimer.h sim.h setup.h

default: $(APP_NAME)-mpi $(APP_NAME)-seq

.PHONY: convert clean

convert:
	mogrify -format jpg out/*.ppm
	rm out/*.ppm

clean:
	rm -rf out/*
	rm $(APP_NAME)-seq
	rm $(APP_NAME)-mpi
	rm -rf sequential

$(APP_NAME)-seq: $(DEPS) $(CXXFILES)
		$(CXX) $(CXXFLAGS) -o $@ $(CXXFILES) -lm

$(APP_NAME)-mpi: $(DEPS) $(CXXFILES)
		$(MPICC) $(CXXFLAGS) $(MPI) -o $@ $(CXXFILES) -lm

