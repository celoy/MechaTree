#$preamble

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = mpif90 #gfortran
# flags for debugging or for maximum performance, comment as necessary
# FCFLAGS = -Ofast -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000
FCFLAGS = -Ofast -fopenmp
# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = EvoluAlgo Forest tree

# "make" builds all
all: $(PROGRAMS)

#$intro
# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

mod_tools.o: precis_mod.o
mod_tree.o: precis_mod.o mod_tools.o
mod_evolu.o: precis_mod.o mod_tree.o
EvoluAlgo.o: precis_mod.o mod_tools.o mod_evolu.o mod_tree.o 
EvoluAlgo: precis_mod.o mod_tools.o mod_evolu.o mod_tree.o
Forest.o: precis_mod.o mod_tools.o mod_evolu.o mod_tree.o 
Forest: precis_mod.o mod_tools.o mod_evolu.o mod_tree.o
tree.o: precis_mod.o mod_tools.o mod_tree.o
tree: precis_mod.o mod_tools.o mod_tree.o

#$conclusion
# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

