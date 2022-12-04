CC = g++
FLAGS = -std=c++17

SOURCEDIR = .
OBJDIR = .
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp,$(OBJDIR)/%,$(SOURCES))

all: $(OBJECTS)

$(OBJECTS): %: %.cpp
	$(CC) $(FLAGS) -o $@ $<

run: ${OBJECTS}
	./${OBJECTS}

.PHONY: clean
clean: ${OBJECTS}
	rm ${OBJECTS}
