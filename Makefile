CC = g++
CXXFLAGS = -std=c++20  
CPPFLAGS = -Wall -O3 -I include


DOXYFILE= include/Doxyfile
SRCS = src/main.cpp  
HDRS = include/matrix.hpp include/function_implementation.hpp include/friend_functions

OBJS = $(SRCS:.cpp=.o) 

TARGET = program

src/%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o $(TARGET)

.PHONY: clean

clean:
	-rm -f $(TARGET) $(OBJS)

doc:
	doxygen $(DOXYFILE)