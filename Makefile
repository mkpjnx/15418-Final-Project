CXX=g++ -m64
CXXFLAGS=-O3 -Wall
DEPS= sim.h
OBJ = main.o sim.o 


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

sequential: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)


convert:
	mogrify -format jpg *.ppm
	rm *.ppm

clean:
	rm -rf *.jpg
	rm -rf *.o
	rm -rf sequential