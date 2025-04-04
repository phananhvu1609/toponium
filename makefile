# Directories for libraries
LHADIR = /home/vuphananh/packages/LHAPDF/lhapdf-653/
GSLDIR = /usr/local/
LTDIR = /home/vuphananh/packages/LoopTools/LoopTools-2.16/x86_64-Linux
LIBDIR = /lib/gcc/x86_64-linux-gnu/11/

# Parameters for compiles. Use -Wfatal-errors to force compile to show only first error. Omit -Wno-unused-variable, -Wno-unused-but-set-variable to warn about unused variables.
FLAG = -O2 -Wall -I$(LHADIR)/include -I$(GSLDIR)/include -I$(LTDIR)/include
LIBS = $(LIBDIR)/libstdc++.a -lm -L$(LHADIR)/lib -lLHAPDF -L$(GSLDIR)/lib -lgsl -lgslcblas -lm -L$(LTDIR)/lib64 -looptools-quad

# Files to compile
SCRS = Utilities.cpp Formulas.cpp Computations.cpp hadronic.cpp
# Compiled files (replace .cpp in $(SCRS) with .o)
OBJS = $(patsubst %.cpp,%.o,$(SCRS))

all: run

# Combine all to executable 'run'
run: $(OBJS)
	$(LTDIR)/bin/f++-quad $(FLAG) -o $@ $^ $(LIBS)

# Compile each .cpp to .o
%.o: %.cpp
	$(LTDIR)/bin/f++-quad -o $@ -c $(FLAG) $< $(LIBS)

clean:
	rm *.o