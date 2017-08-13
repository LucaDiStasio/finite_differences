CXX = g++
CXXFLAGS = -I/usr/lib/x86_64-linux-gnu/
INCLUDES = 
LD = icc
LDFLAGS := -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_program_options -lboost_system

#/home/ubuntu/workspace/boost
#/home/ubuntu/workspace/boost/stage/lib

SOURCES = finitediff.cpp
HEADERS = finitediff.h
OBJ     = main.o finitediff.o

.SUFFIXES: .cpp

all: derivatives

derivatives: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o derivatives $(LDFLAGS)
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	
finitediff.o: finitediff.cpp finitediff.h
	$(CXX) $(CXXFLAGS) -c finitediff.cpp

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*