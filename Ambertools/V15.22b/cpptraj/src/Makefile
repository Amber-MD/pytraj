# Cpptraj stand-alone make file
include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o) $(TIMER)

include ambpdbfiles

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o) MpiRoutines.o

all: cpptraj$(SFX)

install: cpptraj$(SFX) cpp_ambpdb$(SFX) 
	/bin/mv cpptraj$(SFX) $(CPPTRAJBIN)
	/bin/mv cpp_ambpdb $(CPPTRAJBIN)

findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) AmbPDB.cpp > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

cpptraj$(SFX): $(OBJECTS) pub_fft.o $(READLINE)
	$(CXX)  -o cpptraj$(SFX) $(OBJECTS) pub_fft.o $(READLINE) $(LDFLAGS)

cpp_ambpdb: $(AMBPDB_OBJECTS)
	$(CXX) -o cpp_ambpdb $(AMBPDB_OBJECTS) $(LDFLAGS)

cpp_ambpdb.OMP:
	@echo ""
	@echo "INFO: cpp_ambpdb is not OpenMP-enabled. Skipping."
	@echo ""

cpp_ambpdb.MPI:
	@echo ""
	@echo "INFO: cpp_ambpdb is not MPI-enabled. Skipping."
	@echo ""

$(READLINE): $(READLINE_HOME)/Makefile
	cd $(READLINE_HOME) && $(MAKE)

pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

ReadLine.o: ReadLine.cpp
	$(CXX) $(CXXFLAGS) -I$(READLINE_HOME) -c -o $@ ReadLine.cpp

clean:
	/bin/rm -f $(OBJECTS) pub_fft.o cpptraj$(SFX) AmbPDB.o cpp_ambpdb
	cd $(READLINE_HOME) && $(MAKE) clean

uninstall:
	/bin/rm -f $(CPPTRAJBIN)/cpptraj $(CPPTRAJBIN)/cpp_ambpdb readline/libreadline.a

# Header dependencies
include cpptrajdepend
