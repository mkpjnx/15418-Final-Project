CXX=g++ -m64
CXXFLAGS=-Icommon -Iobjs/ -O0 -Wall
OMP=-fopenmp -DOMP

APP_NAME=grayscott
OBJDIR=objs

CXXFILES = main.cpp sim.cpp instrument.cpp
DEPS = cycleTimer.h sim.h

default: $(APP_NAME)-omp $(APP_NAME)-seq

.PHONY: convert clean

convert:
	mogrify -format jpg out/*.ppm
	rm out/*.ppm

clean:
	rm -rf out/*.jpg
	rm -rf out/*.ppm
	rm $(APP_NAME)-seq
	rm $(APP_NAME)-omp
	rm -rf sequential

$(APP_NAME)-seq: $(DEPS) $(CXXFILES)
		$(CXX) $(CXXFLAGS) -o $@ $(CXXFILES) -lm

$(APP_NAME)-omp: $(DEPS) $(CXXFILES)
		$(CXX) $(CXXFLAGS) $(OMP) -o $@ $(CXXFILES) -lm

