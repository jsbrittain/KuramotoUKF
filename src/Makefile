CC = g++
FLAGS = -std=c++17 -g3 -Wall -Wextra -Wno-reorder -c

SOURCEDIR = .
BUILDDIR = ../build
OBJDIR = $(BUILDDIR)/obj

EXECUTABLE = kukf
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

all: dir $(BUILDDIR)/$(EXECUTABLE) tags

dir:
	mkdir -p $(BUILDDIR) $(OBJDIR)

$(BUILDDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) -g $^ -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SOURCEDIR)/%.cpp
	$(CC) $(FLAGS) $< -o $@

tags:
	ctags *.cpp *.hpp

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BUILDDIR)/$(EXECUTABLE) tags
