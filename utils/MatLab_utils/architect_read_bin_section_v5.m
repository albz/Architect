% *** read Architect binary SECTION output ***
% 
% input[1]  -> name with full path
% output[1] -> the entire output
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : A. Marocchino
% Purpose       : read SECTION binary output from Architect: option 5 with Ionisation
% Last modified : 29/8/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kp,wp,dist,n0,dist_um,Nr,Nz,r_mesh,z_mesh,rho_bunch,rho_background,Er,Er_background,Er_bunch,Ez,Ez_background,Ez_bunch,Bphi,Bphi_background,Bphi_bunch,B_ex_poloidal,Jr_bunch,Jr_background,Jz_bunch,Jz_background,Zstar,rho_ions] = architect_read_bin_section_v5(full_name)

file = fopen(full_name,'rb','l');

% --- Output version
output_version=fread(file,1,'int32');

% --- parameters
kp      = fread(file,1,'float64');
wp      = fread(file,1,'float64');
dist    = fread(file,1,'float64');
n0      = fread(file,1,'float64');
dist_um = fread(file,1,'float64');

% - domain axes
Nr=fread(file,1,'int32');                       % half-plane axis dim
Nz=fread(file,1,'int32');                       % z dim
r_mesh =fread(file,Nr,'float64');
z_mesh =fread(file,Nz,'float64');

% - bunch density
rho_bunch=fread(file,Nr*Nz,'float64');              
rho_bunch=reshape(rho_bunch,[Nz Nr]);
rho_bunch=rho_bunch';

% - plasma electron density
rho_background=fread(file,Nr*Nz,'float64');              
rho_background=reshape(rho_background,[Nz Nr]);
rho_background=rho_background';

% - Electric field, transverse
Er=fread(file,Nr*Nz,'float64');
Er=reshape(Er,[Nz Nr]);
Er=Er';
Er_background=fread(file,Nr*Nz,'float64');
Er_background=reshape(Er_background,[Nz Nr]);
Er_background=Er_background';
Er_bunch=fread(file,Nr*Nz,'float64');
Er_bunch=reshape(Er_bunch,[Nz Nr]);
Er_bunch=Er_bunch';

% - Electric field, longitudinal
Ez=fread(file,Nr*Nz,'float64');              
Ez=reshape(Ez,[Nz Nr]);
Ez=Ez';
Ez_background=fread(file,Nr*Nz,'float64');              
Ez_background=reshape(Ez_background,[Nz Nr]);
Ez_background=Ez_background';
Ez_bunch=fread(file,Nr*Nz,'float64');              
Ez_bunch=reshape(Ez_bunch,[Nz Nr]);
Ez_bunch=Ez_bunch';

% - Magnetic field, azimuthal
Bphi=fread(file,Nr*Nz,'float64');              
Bphi=reshape(Bphi,[Nz Nr]);
Bphi=Bphi';
Bphi_background=fread(file,Nr*Nz,'float64');              
Bphi_background=reshape(Bphi_background,[Nz Nr]);
Bphi_background=Bphi_background';
Bphi_bunch=fread(file,Nr*Nz,'float64');              
Bphi_bunch=reshape(Bphi_bunch,[Nz Nr]);
Bphi_bunch=Bphi_bunch';
B_ex_poloidal=fread(file,Nr*Nz,'float64');              
B_ex_poloidal=reshape(B_ex_poloidal,[Nz Nr]);
B_ex_poloidal=B_ex_poloidal';

% - Bunch current density, transverse
Jr_bunch=fread(file,Nr*Nz,'float64');              
Jr_bunch=reshape(Jr_bunch,[Nz Nr]);
Jr_bunch=Jr_bunch';

% - Plasma current density, transverse
Jr_background=fread(file,Nr*Nz,'float64');              
Jr_background=reshape(Jr_background,[Nz Nr]);
Jr_background=Jr_background';

% - Bunch current density, longitudinal
Jz_bunch=fread(file,Nr*Nz,'float64');              
Jz_bunch=reshape(Jz_bunch,[Nz Nr]);
Jz_bunch=Jz_bunch';

% - Plasma current density, longitudinal
Jz_background=fread(file,Nr*Nz,'float64');              
Jz_background=reshape(Jz_background,[Nz Nr]);
Jz_background=Jz_background';

% - Zstar
Zstar=fread(file,Nr*Nz,'float64');              
Zstar=reshape(Zstar,[Nz Nr]);
Zstar=Zstar';

% - rho_i
rho_ions=fread(file,Nr*Nz,'float64');              
rho_ions=reshape(rho_ions,[Nz Nr]);
rho_ions=rho_ions';

fclose('all');
