CC = g++
CXXFLAGS = -std=c++20 -Wall 

DOXYFILE=../../pacs-examples/Examples/DoxyfileCommon
SRCS = src/main.cpp src/matrix.cpp 
HDRS = src/matrix.hpp

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