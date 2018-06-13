CXX = g++
INCLUDEFLAGS = -IH:/GitHub/LATfield2 -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include/x64" -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
LIBFLAGS = -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi
CXXFLAGS = -Wall -g -O2 $(INCLUDEFLAGS) $(LIBFLAGS)
RM = rm -rf
MD = mkdir
TARGET = bin/fieldMath
BUILDDIR = build
SRCDIR = src
SRCEXT = cpp

SOURCES = fieldMath.$(SRCEXT)
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

bin/fieldMath: build/fieldMath.o
	$(CXX) $< -o $@ $(CXXFLAGS)

build/fieldMath.o: src/fieldMath.cpp
	$(CXX) $< -c -o $@ $(CXXFLAGS)

# $(TARGET): $(BUILDDIR)/$(OBJECTS)
# 	$(CXX) $< -o $@ $(CXXFLAGS)

# $(BUILDDIR)/$(OBJECTS): $(SRCDIR)/$(SOURCES)
# 	if not exist $(BUILDDIR) mkdir $(BUILDDIR)
# 	g++ $< -c -o $@ $(CXXFLAGS)

clean:
	$(RM) $(BUILDDIR)
	$(MD) $(BUILDDIR)