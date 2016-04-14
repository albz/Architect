FC          = gfortran-mp-4.7
OPTFC       = -ffree-line-length-none -fmax-stack-var-size=50000000
SRC_FOLDER  = src
OBJ_FOLDER  = .
EXE_FOLDER  = bin
EXE         = Architect

FILES  =  my_types.f90 \
					use_types.f90 \
					Particle_Classes.f90 \
					read_input.f90 \
					fileflags.f90 \
					bunch_generation.f90 \
					utility.f90 \
					moments.f90 \
					bunch_initialization.f90 \
					mesh_generator.f90 \
					linear_algebra_utilities.f90 \
					compute_bunch_current.f90 \
					Fields_init.f90 \
					update_EBfields_leapfrog.f90 \
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


SOURCES     = $(addprefix $(SRC_FOLDER), $(FILES))
OBJECTS     = $(addsuffix .o, $(addprefix $(OBJ_FOLDER)/, $(basename $(FILES))))
MODULES     = $(addsuffix .mod, $(addprefix $(OBJ_FOLDER)/, $(basename $(FILES))))
EXECUTABLE  = $(addprefix $(EXE_FOLDER)/, $(EXE))

all: $(MODULES) $(OBJECTS)
	mkdir -p $(EXE_FOLDER)
	mkdir -p $(OBJ_FOLDER)
	$(FC) $(OPTFC) $(OBJECTS) -o $(EXECUTABLE)




#--- ---#
ifort : FC = ifort
ifort : OPTFC = -O3 -Wl,-stack_size,0x50000000
ifort : all

gflinux : FC = gfortran
gflinux : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gflinux : all




#--- ---#
$(OBJ_FOLDER)/my_types.o: $(SRC_FOLDER)/my_types.F90
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/my_types.mod: $(SRC_FOLDER)/my_types.F90 $(OBJ_FOLDER)/my_types.o
	@true

$(OBJ_FOLDER)/use_types.o: $(SRC_FOLDER)/use_types.F90 $(OBJ_FOLDER)/my_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/use_types.mod: $(SRC_FOLDER)/use_types.F90 $(OBJ_FOLDER)/use_types.o
	@true

$(OBJ_FOLDER)/Particle_Classes.o: $(SRC_FOLDER)/Particle_Classes.F90
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/Particle_Classes.mod: $(SRC_FOLDER)/Particle_Classes.F90 $(OBJ_FOLDER)/Particle_Classes.o
	@true

$(OBJ_FOLDER)/read_input.o: $(SRC_FOLDER)/read_input.F90 \
												    $(OBJ_FOLDER)/my_types.mod \
														$(OBJ_FOLDER)/use_types.mod \
														$(OBJ_FOLDER)/pstruct_data.mod \
														$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/read_input.mod: $(SRC_FOLDER)/read_input.F90 $(OBJ_FOLDER)/read_input.o
	@true

$(OBJ_FOLDER)/fileflags.o: $(SRC_FOLDER)/fileflags.F90
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/fileflags.mod: $(SRC_FOLDER)/fileflags.F90 $(OBJ_FOLDER)/fileflags.o
	@true

$(OBJ_FOLDER)/bunch_generation.o: $(SRC_FOLDER)/bunch_generation.F90
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_generation.mod: $(SRC_FOLDER)/bunch_generation.F90 $(OBJ_FOLDER)/bunch_generation.o
	@true

$(OBJ_FOLDER)/utility.o: $(SRC_FOLDER)/utility.F90 \
												 $(OBJ_FOLDER)/my_types.mod \
												 $(OBJ_FOLDER)/use_types.mod \
												 $(OBJ_FOLDER)/pstruct_data.mod \
												 $(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/utility.mod: $(SRC_FOLDER)/utility.F90 $(OBJ_FOLDER)/utility.o
	@true

$(OBJ_FOLDER)/moments.o: $(SRC_FOLDER)/moments.F90 \
												 $(OBJ_FOLDER)/my_types.mod \
												 $(OBJ_FOLDER)/use_types.mod \
												 $(OBJ_FOLDER)/pstruct_data.mod \
												 $(OBJ_FOLDER)/utilities.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/moments.mod: $(SRC_FOLDER)/moments.F90 $(OBJ_FOLDER)/moments.o
	@true

$(OBJ_FOLDER)/bunch_initialization.o: $(SRC_FOLDER)/bunch_initialization.F90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod \
																			$(OBJ_FOLDER)/pstruct_data.mod \
																			$(OBJ_FOLDER)/architect_class_structure.mod \
																			$(OBJ_FOLDER)/moments.mod \
																			$(OBJ_FOLDER)/bunch_generation.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_initialization.mod: $(SRC_FOLDER)/bunch_initialization.F90 $(OBJ_FOLDER)/bunch_initialization.o
	@true

$(OBJ_FOLDER)/mesh_generator.o: $(SRC_FOLDER)/mesh_generator.F90 \
																$(OBJ_FOLDER)/my_types.mod \
																$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/mesh_generator.mod: $(SRC_FOLDER)/mesh_generator.F90 $(OBJ_FOLDER)/mesh_generator.o
	@true

$(OBJ_FOLDER)/linear_algebra_utilities.o: $(SRC_FOLDER)/linear_algebra_utilities.F90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/linear_algebra_utilities.mod: $(SRC_FOLDER)/linear_algebra_utilities.F90 $(OBJ_FOLDER)/linear_algebra_utilities.o
	@true

$(OBJ_FOLDER)/compute_bunch_current.o: $(SRC_FOLDER)/compute_bunch_current.F90 \
																		   $(OBJ_FOLDER)/my_types.mod \
																			 $(OBJ_FOLDER)/use_types.mod \
																			 $(OBJ_FOLDER)/pstruct_data.mod \
																			 $(OBJ_FOLDER)/architect_class_structure.mod \
																			 $(OBJ_FOLDER)/utilities.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_bunch_current.mod: $(SRC_FOLDER)/compute_bunch_current.F90 $(OBJ_FOLDER)/compute_bunch_current.o
	@true

$(OBJ_FOLDER)/Fields_init.o: $(SRC_FOLDER)/Fields_init.F90 \
														 $(OBJ_FOLDER)/my_types.mod \
														 $(OBJ_FOLDER)/use_types.mod \
														 $(OBJ_FOLDER)/linear_algebra.mod \
														 $(OBJ_FOLDER)/moments.mod \
														 $(OBJ_FOLDER)/Compute_beam_current_FDTD.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/Fields_init.mod: $(SRC_FOLDER)/Fields_init.F90 $(OBJ_FOLDER)/Fields_init.o
	@true

$(OBJ_FOLDER)/update_EBfields_leapfrog.o: $(SRC_FOLDER)/update_EBfields_leapfrog.F90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_EBfields_leapfrog.mod: $(SRC_FOLDER)/update_EBfields_leapfrog.F90 $(OBJ_FOLDER)/update_EBfields_leapfrog.o
	@true

$(OBJ_FOLDER)/update_fluid_upwind.o: $(SRC_FOLDER)/update_fluid_upwind.F90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_upwind.mod: $(SRC_FOLDER)/update_fluid_upwind.F90 $(OBJ_FOLDER)/update_fluid_upwind.o
	@true

$(OBJ_FOLDER)/update_fluid_fct.o: $(SRC_FOLDER)/update_fluid_fct.F90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_fct.mod: $(SRC_FOLDER)/update_fluid_fct.F90 $(OBJ_FOLDER)/update_fluid_fct.o
	@true

$(OBJ_FOLDER)/compute_background_current.o: $(SRC_FOLDER)/compute_background_current.F90 \
																						$(OBJ_FOLDER)/my_types.mod \
																						$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_background_current.mod: $(SRC_FOLDER)/compute_background_current.F90 $(OBJ_FOLDER)/compute_background_current.o
	@true

$(OBJ_FOLDER)/compute_current_manager.o: $(SRC_FOLDER)/compute_current_manager.F90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod \
																					$(OBJ_FOLDER)/utilities.mod \
																					$(OBJ_FOLDER)/pstruct_data.mod \
																					$(OBJ_FOLDER)/architect_class_structure.mod \
																					$(OBJ_FOLDER)/Compute_beam_current_FDTD.mod \
																					$(OBJ_FOLDER)/compute_plasma_current.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_current_manager.mod: $(SRC_FOLDER)/compute_current_manager.F90 $(OBJ_FOLDER)/compute_current_manager.o
	@true

$(OBJ_FOLDER)/update_fluid_manager.o: $(SRC_FOLDER)/update_fluid_manager.F90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_manager.mod: $(SRC_FOLDER)/update_fluid_manager.F90 $(OBJ_FOLDER)/update_fluid_manager.o
	@true

$(OBJ_FOLDER)/particle_pusher.o: $(SRC_FOLDER)/particle_pusher.F90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod \
																	$(OBJ_FOLDER)/utilities.mod \
																	$(OBJ_FOLDER)/pstruct_data.mod \
																	$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/particle_pusher.mod: $(SRC_FOLDER)/particle_pusher.F90 $(OBJ_FOLDER)/particle_pusher.o
	@true

$(OBJ_FOLDER)/window_shifting.o: $(SRC_FOLDER)/window_shifting.F90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod \
																	$(OBJ_FOLDER)/moments.mod \
																	$(OBJ_FOLDER)/pstruct_data.mod \
																	$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/window_shifting.mod: $(SRC_FOLDER)/window_shifting.F90 $(OBJ_FOLDER)/window_shifting.o
	@true

$(OBJ_FOLDER)/new_diagnostics.o: $(SRC_FOLDER)/new_diagnostics.F90 $(OBJ_FOLDER)/my_types.mod $(OBJ_FOLDER)/pstruct_data.mod $(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/new_diagnostics.mod: $(SRC_FOLDER)/new_diagnostics.F90 $(OBJ_FOLDER)/new_diagnostics.o
	@true

$(OBJ_FOLDER)/data_dump.o: $(SRC_FOLDER)/data_dump.F90 $(OBJ_FOLDER)/my_types.mod $(OBJ_FOLDER)/use_types.mod $(OBJ_FOLDER)/pstruct_data.mod $(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/data_dump.mod: $(SRC_FOLDER)/data_dump.F90 $(OBJ_FOLDER)/data_dump.o
	@true

$(OBJ_FOLDER)/architect.o: $(SRC_FOLDER)/architect.F90 \
													 $(OBJ_FOLDER)/my_types.mod \
													 $(OBJ_FOLDER)/use_types.mod \
													 $(OBJ_FOLDER)/read_input_module.mod \
													 $(OBJ_FOLDER)/make_a_mesh.mod \
													 $(OBJ_FOLDER)/diagnostics_on_bunches.mod \
													 $(OBJ_FOLDER)/fields_FDTD.mod \
													 $(OBJ_FOLDER)/move_window_FDTD.mod \
													 $(OBJ_FOLDER)/computecurrentFDTD.mod \
													 $(OBJ_FOLDER)/Compute_beam_current_FDTD.mod \
													 $(OBJ_FOLDER)/moveparticle_FDTD.mod \
													 $(OBJ_FOLDER)/fluidadvance_FDTD.mod \
													 $(OBJ_FOLDER)/bunch_generation.mod \
													 $(OBJ_FOLDER)/data_dumping.mod \
													 $(OBJ_FOLDER)/pstruct_data.mod \
													 $(OBJ_FOLDER)/architect_class_structure.mod \
													 $(OBJ_FOLDER)/moments.mod \
													 $(OBJ_FOLDER)/initialize_bunch.mod \
													 $(OBJ_FOLDER)/utilities.mod \
													 $(OBJ_FOLDER)/linear_algebra.mod \
													 $(OBJ_FOLDER)/init_fields.mod
	$(FC) $(OPTFC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/architect.mod: $(SRC_FOLDER)/architect.F90 $(OBJ_FOLDER)/architect.o
	@true


clean:
	@echo 'Cleaning all object files'
	rm -f $(OBJECTS)
	rm -f *.mod

cleanall:
	@echo 'Cleaning all object files, subdirs and EXECutable'
	rm -f $(OBJECTS) $(MODULES) *.mod $(EXECUTABLE) opt_report.txt *.exe *~ .*~ log.txt *.optrpt bin/*.optrpt
	rm -rf $(EXE_FOLDER)
