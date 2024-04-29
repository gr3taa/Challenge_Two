CC = g++
CFLAGS = -std=c++20 -Wall -g
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Lista dei file oggetto
OBJECTS = $(OBJ_DIR)/main.o $(OBJ_DIR)/matrix.o

# Nome del programma finale
TARGET = $(BIN_DIR)/myprogram

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS)

# Compila il file main.cpp
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/helper.hpp $(SRC_DIR)/matrix.hpp $(SRC_DIR)/matrix.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Compila il file matrix.cpp
$(OBJ_DIR)/matrix.o: $(SRC_DIR)/matrix.cpp $(SRC_DIR)/matrix.hpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean
