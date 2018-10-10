# Guide for **particle data structure** <br/> for 3D and 2D-cylindrical bunches

## data structure (3D Cartesian Bunch):
Particles coordinates and IDs are stored in the *14 components* of the bunch class. The bunch class is so identified **bunch(bunch_number)%part(particle_number)%cmp(component_number)**. The 14 components, or cmp, have the following meaning:  
+ cmp(1) >> *x* coordinate (in µm)
+ cmp(2) >> *y* coordinate (in µm)
+ cmp(3) >> *z* coordinate (in µm)
+ cmp(4) >> *px* momentum (in beta*gamma)
+ cmp(5) >> *py* momentum (in beta*gamma)
+ cmp(6) >> *pz* momentum (in beta*gamma)


+ cmp(7) >> *cut* component. When cmp(7)=0 the cut is activates and the code does not track the particle anymore
+ cmp(8) >> *diagnostic cut* (*dcut*). Similar to cmp(7) but in this case


+ cmp(9)  >> *x_old or x(t-1)* past coordinate (in µm)
+ cmp(10)  >> *y_old or y(t-1)* past coordinate (in µm)
+ cmp(11)  >> *z_old or z(t-1)* past coordinate (in µm)


+ cmp(12)  >> *particle in charge* (in nC)
+ cmp(13)  >> *number of electrons for macroparticle)*


+ cmp(14)  >> *diagnostic flag* this flag has been added for diagnostic and it is constantly changed within the code. If the value of cmp(14) is 0 it is not counted in any integrated diagnostic, if it is 1, on the contrary, it is considered within the integrated diagnostic

## Re-Usage of the data structure for the 2D Cylindrical symmetric Bunch case:
