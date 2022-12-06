CC = g++
FLAGS = -std=c++17 -g -c

SOURCEDIR = .
BUILDDIR = ../build
OBJDIR = $(BUILDDIR)/obj

EXECUTABLE = kukf
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

all: dir $(BUILDDIR)/$(EXECUTABLE)

dir:
	mkdir -p $(BUILDDIR) $(OBJDIR)

debug: FLAGS += -DDEBUG -g

$(BUILDDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) -g $^ -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SOURCEDIR)/%.cpp
	$(CC) $(FLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BUILDDIR)/$(EXECUTABLE)
