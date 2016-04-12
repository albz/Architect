FF       = ifort
OPT      = -O3 -Wl,-stack_size,0x50000000
EXEC     = Architect

SOURCES  = my_types.f90 use_types.f90 Particle_Classes.f90 read_input.f90 \
        fileflags.f90 bunch_generation.f90 utility.f90 moments.f90 \
	bunch_initialization.f90 \
	mesh_generator.f90 \
	linear_algebra_utilities.f90 \
	compute_bunch_current.f90 \
	Fields_init.f90  update_EBfields_leapfrog.f90 \
	update_fluid_upwind.f90 \
	update_fluid_fct.f90 \
	compute_background_current.f90 \
	compute_current_manager.f90 \
	update_fluid_manager.f90 \
	particle_pusher.f90 \
	window_shifting.f90 \
	new_diagnostics.f90 \
	data_dump.f90 \
	architect.f90


OBJECTS  = $(addsuffix .o, $(basename $(SOURCES)))

all       : $(EXEC)

$(EXEC)   : $(OBJECTS)
	$(FF) $(OPT) -o $(EXEC) $(OBJECTS)


gfortran47 : FF = gfortran-mp-4.7
gfortran47 : OPT = -ffree-line-length-none -fmax-stack-var-size=50000000
gfortran47 : $(EXEC)

gfortran48 : FF = gfortran-mp-4.5
gfortran48 : OPT = -ffree-line-length-none #-fmax-stack-var-size=50000000
gfortran48 : $(EXEC)

gfortran5 : FF = gfortran-5
gfortran5 : OPT = -ffree-line-length-none #-fmax-stack-var-size=50000000
gfortran5 : $(EXEC)

gfortran-lisa: FF = gfortran
gfortran-lisa: OPT = -O3 -ffree-line-length-none -fopenmp -cpp
gfortran-lisa: $(EXEC)

gflinux : FF = gfortran
gflinux : OPT = -ffree-line-length-none -fmax-stack-var-size=50000000
gflinux : $(EXEC)

clean :
	echo 'Cleaning all object files'
	rm -f *.mod *.o

#compilation rule(s)
%.o : %.f90
	$(FF) $(OPT) -c $<
