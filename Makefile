
CXX = g++
FFTWDIR = /sw/fftw/3.3.4-gcc8
CXXFLAGSINIT = -std=c++17 -O3 -malign-double   
CXXFLAGSDYN = -std=c++17 -O3 -malign-double -I$(FFTWDIR)/include # -fopenmp 
GCCOMPFLAGS = -fopenmp -liomp5 -lirc -lpthread
CXXFLAGS = -std=c++17 -O3 -malign-double -fopenmp 

#############################

INCDIR := ./include
SRCDIR := ./src
TESTDIR := ./test
OBJDIR := ./obj
BINDIR := ./bin
OUTPUTDIR := ./output_data
LDFLAGS := -L$(FFTWDIR)/lib -lfftw3
LDFLAGSPARALLEL := -L$(FFTWDIR)/lib -lfftw3_omp -lfftw3

#############################

OBJDIRINIT := $(OBJDIR)/objinit
OBJDIRDYN := $(OBJDIR)/objdyn
OBJDIRTEST := $(OBJDIR)/objtest

############################

# Header files for different executions

INITHEADER = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp $(INCDIR)/SOFT/tunneling_init_state.hpp  
DYNHEADER  = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp $(INCDIR)/SOFT/TROTTER_prop.hpp $(INCDIR)/SOFT/data_process.hpp
TESTHEADER = $(INCDIR)/SOFT/constants.hpp $(INCDIR)/SOFT/fileops.hpp

#############################

# Initial state preparation

init.exe: $(OBJDIRINIT)/constants.o $(OBJDIRINIT)/fileops.o $(OBJDIRINIT)/tunneling_init_state.o $(OBJDIRINIT)/get_init_state.o
	@echo "Generating Initial State ..."
	$(CXX) -o $@	$^ $(GCCOMPFLAGS) 

$(OBJDIRINIT)/constants.o: $(SRCDIR)/constants.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGSINIT) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/fileops.o: $(SRCDIR)/fileops.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGSINIT) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/tunneling_init_state.o: $(SRCDIR)/tunneling_init_state.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGSINIT) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRINIT)/get_init_state.o: $(SRCDIR)/get_init_state.cpp $(INITHEADER)
	$(CXX) $(CXXFLAGSINIT) -I$(INCDIR) $(GCCOMPFLAGS) -c $< -o $@ 


#############################

# Running SOFT dynamics

dyn.exe: $(OBJDIRDYN)/constants.o $(OBJDIRDYN)/fileops.o $(OBJDIRDYN)/TROTTER_prop.o $(OBJDIRDYN)/dynamics.o $(OBJDIRDYN)/data_process.o
	@echo "Running Split-Operator Fourier Transform dynamics ..."
	$(CXX) -o $@	$^	$(LDFLAGSPARALLEL) $(GCCOMPFLAGS)

$(OBJDIRDYN)/constants.o: $(SRCDIR)/constants.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGSDYN) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/fileops.o: $(SRCDIR)/fileops.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGSDYN) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/TROTTER_prop.o: $(SRCDIR)/TROTTER_prop.cpp $(DYNHEADER) 
	$(CXX) $(CXXFLAGSDYN) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRDYN)/data_process.o: $(SRCDIR)/data_process.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGSDYN) -I$(INCDIR) $(GCCOMPFLAGS) -c $< -o $@ 

$(OBJDIRDYN)/dynamics.o: $(SRCDIR)/dynamics.cpp $(DYNHEADER)
	$(CXX) $(CXXFLAGSDYN) -I$(INCDIR) $(GCCOMPFLAGS) -c $< -o $@ 


#############################

# Running analytical dynamics

testanal.exe: $(OBJDIRTEST)/constants.o $(OBJDIRTEST)/fileops.o $(OBJDIRTEST)/dynamics_analytical.o	
	@echo "Testing .."
	@echo "Running analytical time-evolution for initial state from linear combination of eigenstates"
	$(CXX)	-o $@	$^


$(OBJDIRTEST)/constants.o: $(TESTDIR)/constants.cpp $(TESTHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRTEST)/fileops.o: $(TESTDIR)/fileops.cpp $(TESTHEADER)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

$(OBJDIRTEST)/dynamics_analytical.o: $(TESTDIR)/dynamics_analytical.cpp 
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ 

.PHONY: clean


move:
	@echo "Moving execs after a run to bin directory ..."
	mv -f *.exe ./bin


clean:
	@echo "Cleaning all exec, objects, binaries in output directory ..."
	rm -f $(OBJDIRINIT)/* 
	rm -f $(OBJDIRDYN)/* 
	rm -f $(OBJDIRTEST)/* 
	rm -f $(BINDIR)/*
	rm -f $(OUTPUTDIR)/*.bin
	
