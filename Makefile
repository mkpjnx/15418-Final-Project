CXX=g++ -m64
CXXFLAGS=-Icommon -Iobjs/ -O3 -Wall

APP_NAME=grayscott
OBJDIR=objs

OBJ = main.o sim.o instrument.o

default: $(APP_NAME)

.PHONY: dirs convert clean

dirs:
		/bin/mkdir -p $(OBJDIR)/

convert:
	mogrify -format jpg out/*.ppm
	rm out/*.ppm

clean:
	rm -rf out/*.jpg
	rm -rf out/*.ppm
	rm $(APP_NAME)
	rm -rf $(OBJDIR)/*.o
	rm -rf sequential

OBJS = $(OBJDIR)/main.o $(OBJDIR)/sim.o $(OBJDIR)/instrument.o

$(APP_NAME): dirs $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm

$(OBJDIR)/%.o: %.cpp
		$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/main.o: cycleTimer.h
$(OBJDIR)/instrument.o: cycleTimer.h