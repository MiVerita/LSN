CXX       = mpic++
CXXFLAGS  = -std=c++11 -O2
SRCS      = random.cpp System.cpp main.cpp
OBJS      = $(SRCS:.cpp=.o)
EXEC      = main.exe

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS)

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

System.o: System.cpp System.h random.h
	$(CXX) $(CXXFLAGS) -c System.cpp

main.o: main.cpp System.h random.h
	$(CXX) $(CXXFLAGS) -c main.cpp

	