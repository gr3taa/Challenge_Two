CC = g++
CXXFLAGS = -std=c++20  
CPPFLAGS = -Wall -O3 -I include


DOXYFILE= Doxyfile
SRCS = src/main.cpp  
HDRS = src/matrix.hpp src/function_implementation.hpp

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