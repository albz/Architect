FC          = gfortran
CC         = g++
OPTFC       = -ffree-line-length-none -fmax-stack-var-size=50000000
OPTCC       = -O3
SRC_FOLDER  = src
OBJ_FOLDER  = obj
EXE_FOLDER  = bin
EXE         = Architect
BOOST_LIB   = /usr/lib/
BOOST_INC   = /usr/include/
BOOST_FS    = -lboost_filesystem
BOOST_S     = -lboost_system
STDCPP_LINK = -lstdc++
MODULE_REDIRECT = -I$(OBJ_FOLDER) -J$(OBJ_FOLDER)

FILES  =  my_types.f90 \
					use_types.f90 \
					random_numbers_functions.f90 \
					shapiro_wilks.f90 \
					Particle_Classes.f90 \
					read_input.f90 \
					fileflags.f90 \
					bunch_generation.f90 \
					utility.f90 \
					bunch_moments.f90 \
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
					ion_background.f90 \
					ionisation.f90 \
					window_shifting.f90 \
					bunch_diagnostics.f90 \
					grid_diagnostics.f90 \
					data_dump.f90 \
					data_dump_xlm.cpp \
					dump_status.f90 \
					architect.f90


SOURCES     = $(addprefix $(SRC_FOLDER), $(FILES))
OBJECTS     = $(addsuffix .o, $(addprefix $(OBJ_FOLDER)/, $(basename $(FILES))))
MODULES     = $(addsuffix .mod, $(addprefix $(OBJ_FOLDER)/, $(basename $(FILES))))
EXECUTABLE  = $(addprefix $(EXE_FOLDER)/, $(EXE))

all: dirtree $(MODULES) $(OBJECTS)
	$(FC) $(OPTFC) -J$(OBJ_FOLDER) $(OBJECTS) -o $(EXECUTABLE)



dirtree:
	mkdir -p $(EXE_FOLDER)
	mkdir -p $(OBJ_FOLDER)



#--- ---#
ifort : FC = ifort
ifort : OPTFC = -O3 -Wl,-stack_size,0x50000000
ifort : MODULE_REDIRECT = -I$(OBJ_FOLDER) -module $(OBJ_FOLDER)
ifort : all

gflinux : FC = gfortran
gflinux : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gflinux : all

gflinux_lxp : FC = gfortran44
gflinux_lxp : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gflinux_lxp : all

gfmp45  : FC = gfortran-mp-4.5
gfmp45  : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gfmp45  : all
gfmp45  : CC = gcc-mp-4.5

gfortran48 : FC = gfortran-mp-4.8
gfortran48 : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gfortran48 : all

gfortran49 : FC = gfortran-4.9
gfortran49 : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gfortran49 : all

gfortran5 : FC = gfortran-5
gfortran5 : OPTFC = -ffree-line-length-none -fmax-stack-var-size=50000000
gfortran5 : all



#--- ---#
$(OBJ_FOLDER)/my_types.o: $(SRC_FOLDER)/my_types.f90
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/my_types.mod: $(SRC_FOLDER)/my_types.f90 $(OBJ_FOLDER)/my_types.o
	@true

$(OBJ_FOLDER)/use_types.o: $(SRC_FOLDER)/use_types.f90 $(OBJ_FOLDER)/my_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/use_types.mod: $(SRC_FOLDER)/use_types.f90 $(OBJ_FOLDER)/use_types.o
	@true

$(OBJ_FOLDER)/random_numbers_functions.o: $(SRC_FOLDER)/random_numbers_functions.f90
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/random_numbers_functions.mod: $(SRC_FOLDER)/random_numbers_functions.f90 $(OBJ_FOLDER)/random_numbers_functions.o
	@true

$(OBJ_FOLDER)/shapiro_wilks.o: $(SRC_FOLDER)/shapiro_wilks.f90
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/shapiro_wilks.mod: $(SRC_FOLDER)/shapiro_wilks.f90 $(OBJ_FOLDER)/shapiro_wilks.o
	@true

$(OBJ_FOLDER)/Particle_Classes.o: $(SRC_FOLDER)/Particle_Classes.f90
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/Particle_Classes.mod: $(SRC_FOLDER)/Particle_Classes.f90 $(OBJ_FOLDER)/Particle_Classes.o
	@true

$(OBJ_FOLDER)/read_input.o: $(SRC_FOLDER)/read_input.f90 \
												    $(OBJ_FOLDER)/my_types.mod \
														$(OBJ_FOLDER)/use_types.mod \
														$(OBJ_FOLDER)/pstruct_data.mod \
														$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/read_input.mod: $(SRC_FOLDER)/read_input.f90 $(OBJ_FOLDER)/read_input.o
	@true

$(OBJ_FOLDER)/fileflags.o: $(SRC_FOLDER)/fileflags.f90
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/fileflags.mod: $(SRC_FOLDER)/fileflags.f90 $(OBJ_FOLDER)/fileflags.o
	@true

$(OBJ_FOLDER)/bunch_generation.o: $(SRC_FOLDER)/bunch_generation.f90 \
																	$(OBJ_FOLDER)/random_numbers_functions.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_generation.mod: $(SRC_FOLDER)/bunch_generation.f90 $(OBJ_FOLDER)/bunch_generation.o
	@true

$(OBJ_FOLDER)/utility.o: $(SRC_FOLDER)/utility.f90 \
												 $(OBJ_FOLDER)/my_types.mod \
												 $(OBJ_FOLDER)/use_types.mod \
												 $(OBJ_FOLDER)/pstruct_data.mod \
												 $(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/utility.mod: $(SRC_FOLDER)/utility.f90 $(OBJ_FOLDER)/utility.o
	@true

$(OBJ_FOLDER)/bunch_moments.o: $(SRC_FOLDER)/bunch_moments.f90 \
												 $(OBJ_FOLDER)/my_types.mod \
												 $(OBJ_FOLDER)/use_types.mod \
												 $(OBJ_FOLDER)/pstruct_data.mod \
												 $(OBJ_FOLDER)/utilities.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_moments.mod: $(SRC_FOLDER)/bunch_moments.f90 $(OBJ_FOLDER)/bunch_moments.o
	@true

$(OBJ_FOLDER)/bunch_initialization.o: $(SRC_FOLDER)/bunch_initialization.f90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod \
																			$(OBJ_FOLDER)/pstruct_data.mod \
																			$(OBJ_FOLDER)/architect_class_structure.mod \
																			$(OBJ_FOLDER)/bunch_moments.mod \
																			$(OBJ_FOLDER)/bunch_generation.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_initialization.mod: $(SRC_FOLDER)/bunch_initialization.f90 $(OBJ_FOLDER)/bunch_initialization.o
	@true

$(OBJ_FOLDER)/mesh_generator.o: $(SRC_FOLDER)/mesh_generator.f90 \
																$(OBJ_FOLDER)/my_types.mod \
																$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/mesh_generator.mod: $(SRC_FOLDER)/mesh_generator.f90 $(OBJ_FOLDER)/mesh_generator.o
	@true

$(OBJ_FOLDER)/linear_algebra_utilities.o: $(SRC_FOLDER)/linear_algebra_utilities.f90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/linear_algebra_utilities.mod: $(SRC_FOLDER)/linear_algebra_utilities.f90 $(OBJ_FOLDER)/linear_algebra_utilities.o
	@true

$(OBJ_FOLDER)/compute_bunch_current.o: $(SRC_FOLDER)/compute_bunch_current.f90 \
																		   $(OBJ_FOLDER)/my_types.mod \
																			 $(OBJ_FOLDER)/use_types.mod \
																			 $(OBJ_FOLDER)/pstruct_data.mod \
																			 $(OBJ_FOLDER)/architect_class_structure.mod \
																			 $(OBJ_FOLDER)/utilities.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_bunch_current.mod: $(SRC_FOLDER)/compute_bunch_current.f90 $(OBJ_FOLDER)/compute_bunch_current.o
	@true

$(OBJ_FOLDER)/Fields_init.o: $(SRC_FOLDER)/Fields_init.f90 \
														 $(OBJ_FOLDER)/my_types.mod \
														 $(OBJ_FOLDER)/use_types.mod \
														 $(OBJ_FOLDER)/linear_algebra.mod \
														 $(OBJ_FOLDER)/bunch_moments.mod \
														 $(OBJ_FOLDER)/compute_beam_current_fdtd.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/Fields_init.mod: $(SRC_FOLDER)/Fields_init.f90 $(OBJ_FOLDER)/Fields_init.o
	@true

$(OBJ_FOLDER)/update_EBfields_leapfrog.o: $(SRC_FOLDER)/update_EBfields_leapfrog.f90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_EBfields_leapfrog.mod: $(SRC_FOLDER)/update_EBfields_leapfrog.f90 $(OBJ_FOLDER)/update_EBfields_leapfrog.o
	@true

$(OBJ_FOLDER)/update_fluid_upwind.o: $(SRC_FOLDER)/update_fluid_upwind.f90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_upwind.mod: $(SRC_FOLDER)/update_fluid_upwind.f90 $(OBJ_FOLDER)/update_fluid_upwind.o
	@true

$(OBJ_FOLDER)/update_fluid_fct.o: $(SRC_FOLDER)/update_fluid_fct.f90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_fct.mod: $(SRC_FOLDER)/update_fluid_fct.f90 $(OBJ_FOLDER)/update_fluid_fct.o
	@true

$(OBJ_FOLDER)/compute_background_current.o: $(SRC_FOLDER)/compute_background_current.f90 \
																						$(OBJ_FOLDER)/my_types.mod \
																						$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_background_current.mod: $(SRC_FOLDER)/compute_background_current.f90 $(OBJ_FOLDER)/compute_background_current.o
	@true

$(OBJ_FOLDER)/compute_current_manager.o: $(SRC_FOLDER)/compute_current_manager.f90 \
																					$(OBJ_FOLDER)/my_types.mod \
																					$(OBJ_FOLDER)/use_types.mod \
																					$(OBJ_FOLDER)/utilities.mod \
																					$(OBJ_FOLDER)/pstruct_data.mod \
																					$(OBJ_FOLDER)/architect_class_structure.mod \
																					$(OBJ_FOLDER)/compute_beam_current_fdtd.mod \
																					$(OBJ_FOLDER)/compute_plasma_current.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/compute_current_manager.mod: $(SRC_FOLDER)/compute_current_manager.f90 $(OBJ_FOLDER)/compute_current_manager.o
	@true

$(OBJ_FOLDER)/update_fluid_manager.o: $(SRC_FOLDER)/update_fluid_manager.f90 \
																			$(OBJ_FOLDER)/my_types.mod \
																			$(OBJ_FOLDER)/use_types.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/update_fluid_manager.mod: $(SRC_FOLDER)/update_fluid_manager.f90 $(OBJ_FOLDER)/update_fluid_manager.o
	@true

$(OBJ_FOLDER)/particle_pusher.o: $(SRC_FOLDER)/particle_pusher.f90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod \
																	$(OBJ_FOLDER)/utilities.mod \
																	$(OBJ_FOLDER)/pstruct_data.mod \
																	$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/particle_pusher.mod: $(SRC_FOLDER)/particle_pusher.f90 $(OBJ_FOLDER)/particle_pusher.o
	@true

$(OBJ_FOLDER)/ion_background.o: $(SRC_FOLDER)/ion_background.f90 \
																$(OBJ_FOLDER)/my_types.mod \
																$(OBJ_FOLDER)/use_types.mod \
																$(OBJ_FOLDER)/utilities.mod \
																$(OBJ_FOLDER)/pstruct_data.mod \
																$(OBJ_FOLDER)/compute_plasma_current.mod \
																$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/ion_background.mod: $(SRC_FOLDER)/ion_background.f90 $(OBJ_FOLDER)/ion_background.o
		@true

$(OBJ_FOLDER)/ionisation.o: $(SRC_FOLDER)/ionisation.f90 \
														$(OBJ_FOLDER)/my_types.mod \
														$(OBJ_FOLDER)/use_types.mod \
														$(OBJ_FOLDER)/utilities.mod \
														$(OBJ_FOLDER)/pstruct_data.mod \
														$(OBJ_FOLDER)/architect_class_structure.mod \
														$(OBJ_FOLDER)/ion_background.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/ionisation.mod: $(SRC_FOLDER)/ionisation.f90 $(OBJ_FOLDER)/ionisation.o
		@true

$(OBJ_FOLDER)/window_shifting.o: $(SRC_FOLDER)/window_shifting.f90 \
																	$(OBJ_FOLDER)/my_types.mod \
																	$(OBJ_FOLDER)/use_types.mod \
																	$(OBJ_FOLDER)/bunch_moments.mod \
																	$(OBJ_FOLDER)/pstruct_data.mod \
																	$(OBJ_FOLDER)/architect_class_structure.mod \
																	$(OBJ_FOLDER)/compute_plasma_current.mod \
																	$(OBJ_FOLDER)/ion_background.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/window_shifting.mod: $(SRC_FOLDER)/window_shifting.f90 $(OBJ_FOLDER)/window_shifting.o
	@true

$(OBJ_FOLDER)/bunch_diagnostics.o: $(SRC_FOLDER)/bunch_diagnostics.f90 \
																	 $(OBJ_FOLDER)/my_types.mod \
																	 $(OBJ_FOLDER)/pstruct_data.mod \
																	 $(OBJ_FOLDER)/architect_class_structure.mod \
																	 $(OBJ_FOLDER)/random_numbers_functions.mod \
																	 $(OBJ_FOLDER)/shapiro_wilks.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/bunch_diagnostics.mod: $(SRC_FOLDER)/bunch_diagnostics.f90 $(OBJ_FOLDER)/bunch_diagnostics.o
	@true

$(OBJ_FOLDER)/grid_diagnostics.o: $(SRC_FOLDER)/grid_diagnostics.f90 \
																	 $(OBJ_FOLDER)/my_types.mod \
 																	 $(OBJ_FOLDER)/pstruct_data.mod \
																	 $(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/grid_diagnostics.mod: $(SRC_FOLDER)/grid_diagnostics.f90 $(OBJ_FOLDER)/grid_diagnostics.o
	@true

$(OBJ_FOLDER)/data_dump.o: $(SRC_FOLDER)/data_dump.f90 \
	                         $(OBJ_FOLDER)/my_types.mod \
													 $(OBJ_FOLDER)/use_types.mod \
													 $(OBJ_FOLDER)/pstruct_data.mod \
													 $(OBJ_FOLDER)/architect_class_structure.mod \
													 $(OBJ_FOLDER)/ion_background.mod \
													 $(OBJ_FOLDER)/grid_diagnostics.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/data_dump.mod: $(SRC_FOLDER)/data_dump.f90 $(OBJ_FOLDER)/data_dump.o
	@true

$(OBJ_FOLDER)/data_dump_xlm.o: $(SRC_FOLDER)/data_dump_xlm.cpp
	$(CC) $(OPTCC) -I$(BOOST_INC) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/data_dump_xlm.mod: $(SRC_FOLDER)/data_dump_xlm.cpp $(OBJ_FOLDER)/data_dump_xlm.o
	@true

$(OBJ_FOLDER)/dump_status.o: $(SRC_FOLDER)/dump_status.f90 \
															$(OBJ_FOLDER)/my_types.mod \
															$(OBJ_FOLDER)/use_my_types.mod \
															$(OBJ_FOLDER)/pstruct_data.mod \
															$(OBJ_FOLDER)/architect_class_structure.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/dump_status.mod: $(SRC_FOLDER)/dump_status.f90 $(OBJ_FOLDER)/dump_status.o
	@true

$(OBJ_FOLDER)/architect.o: $(SRC_FOLDER)/architect.f90 \
													 $(OBJ_FOLDER)/my_types.mod \
													 $(OBJ_FOLDER)/use_types.mod \
													 $(OBJ_FOLDER)/read_input_module.mod \
													 $(OBJ_FOLDER)/make_a_mesh.mod \
													 $(OBJ_FOLDER)/fields_fdtd.mod \
													 $(OBJ_FOLDER)/move_window_fdtd.mod \
													 $(OBJ_FOLDER)/computecurrentfdtd.mod \
													 $(OBJ_FOLDER)/compute_beam_current_fdtd.mod \
													 $(OBJ_FOLDER)/moveparticle_fdtd.mod \
													 $(OBJ_FOLDER)/fluidadvance_fdtd.mod \
													 $(OBJ_FOLDER)/bunch_generation.mod \
													 $(OBJ_FOLDER)/data_dumping.mod \
													 $(OBJ_FOLDER)/pstruct_data.mod \
													 $(OBJ_FOLDER)/architect_class_structure.mod \
													 $(OBJ_FOLDER)/bunch_moments.mod \
													 $(OBJ_FOLDER)/initialize_bunch.mod \
													 $(OBJ_FOLDER)/utilities.mod \
													 $(OBJ_FOLDER)/linear_algebra.mod \
													 $(OBJ_FOLDER)/init_fields.mod \
													 $(OBJ_FOLDER)/ion_background.mod \
													 $(OBJ_FOLDER)/dump_status.mod \
													 $(OBJ_FOLDER)/ionisation_module.mod
	$(FC) $(OPTFC) $(MODULE_REDIRECT) -c -o $@ $< $(REDIRECT)
$(OBJ_FOLDER)/architect.mod: $(SRC_FOLDER)/architect.f90 $(OBJ_FOLDER)/architect.o
	@true


clean:
	@echo 'Cleaning all object files'
	rm -f $(OBJECTS)
	rm -f *.mod

cleanall:
	@echo 'Cleaning all object files, subdirs and EXECutable'
	rm -f $(OBJECTS) $(MODULES) *.mod $(EXECUTABLE) opt_report.txt *.exe *~ .*~ log.txt *.optrpt bin/*.optrpt
	rm -rf $(EXE_FOLDER) $(OBJ_FOLDER)
