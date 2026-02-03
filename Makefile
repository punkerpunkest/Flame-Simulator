CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2 -Iinclude -I/opt/homebrew/include -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
DEBUGFLAGS = -std=c++20 -Wall -Wextra -g -O0 -Iinclude -I/opt/homebrew/include -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib -lomp
SFML_LIBS = -lsfml-graphics -lsfml-window -lsfml-system

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

TARGET = $(BUILD_DIR)/flamesim
DEBUG_TARGET = $(BUILD_DIR)/flamesim_debug

OBJECTS = $(BUILD_DIR)/main.o $(BUILD_DIR)/FlameSim.o
DEBUG_OBJECTS = $(BUILD_DIR)/main_debug.o $(BUILD_DIR)/FlameSim_debug.o

HEADERS = $(INCLUDE_DIR)/FlameSim.hpp

.PHONY: all debug run run-debug clean rebuild

all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(SFML_LIBS)

debug: $(DEBUG_TARGET)

$(DEBUG_TARGET): $(DEBUG_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(DEBUGFLAGS) $(DEBUG_OBJECTS) -o $(DEBUG_TARGET) $(LDFLAGS) $(SFML_LIBS)

$(BUILD_DIR)/main.o: main.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/FlameSim.o: $(SRC_DIR)/FlameSim.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/main_debug.o: main.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(DEBUGFLAGS) -c $< -o $@

$(BUILD_DIR)/FlameSim_debug.o: $(SRC_DIR)/FlameSim.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(DEBUGFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET)

run-debug: $(DEBUG_TARGET)
	./$(DEBUG_TARGET)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all