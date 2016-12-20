% *** read Architect binary SECTION output ***
% 
% input[1]  -> name with full path
% output[1] -> the entire output
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : A. Marocchino
% Purpose       : read SECTION binary output from Architect: option 4 with Ionisation
% Last modified : 29/8/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist,Nr,Nz,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz,Zstar,rho_i] = architect_read_bin_section_v4(full_name)

file = fopen(full_name,'rb','l');

% --- Output version
output_version=fread(file,1,'int32');

% --- Traveled distance
dist=fread(file,1,'int32');

% - domain axes
Nr=fread(file,1,'int32');                       % half-plane axis dim
Nz=fread(file,1,'int32');                       % z dim
r_mesh =fread(file,Nr,'float64');
z_mesh =fread(file,Nz,'float64');

% - bunch density
rho_b=fread(file,Nr*Nz,'float64');              
rho_b=reshape(rho_b,[Nz Nr]);
rho_b=rho_b';

% - plasma electron density
n_bck=fread(file,Nr*Nz,'float64');              
n_bck=reshape(n_bck,[Nz Nr]);
n_bck=n_bck';

% - Electric field, transverse
Er=fread(file,Nr*Nz,'float64');
Er=reshape(Er,[Nz Nr]);
Er=Er';
Er_bck=fread(file,Nr*Nz,'float64');
Er_bck=reshape(Er_bck,[Nz Nr]);
Er_bck=Er_bck';
Er_b=fread(file,Nr*Nz,'float64');
Er_b=reshape(Er_b,[Nz Nr]);
Er_b=Er_b';

% - Electric field, longitudinal
Ez=fread(file,Nr*Nz,'float64');              
Ez=reshape(Ez,[Nz Nr]);
Ez=Ez';
Ez_bck=fread(file,Nr*Nz,'float64');              
Ez_bck=reshape(Ez_bck,[Nz Nr]);
Ez_bck=Ez_bck';
Ez_b=fread(file,Nr*Nz,'float64');              
Ez_b=reshape(Ez_b,[Nz Nr]);
Ez_b=Ez_b';

% - Magnetic field, azimuthal
Bphi=fread(file,Nr*Nz,'float64');              
Bphi=reshape(Bphi,[Nz Nr]);
Bphi=Bphi';
Bphi_bck=fread(file,Nr*Nz,'float64');              
Bphi_bck=reshape(Bphi_bck,[Nz Nr]);
Bphi_bck=Bphi_bck';
Bphi_b=fread(file,Nr*Nz,'float64');              
Bphi_b=reshape(Bphi_b,[Nz Nr]);
Bphi_b=Bphi_b';
B_ex_poloidal=fread(file,Nr*Nz,'float64');              
B_ex_poloidal=reshape(B_ex_poloidal,[Nz Nr]);
B_ex_poloidal=B_ex_poloidal';

% - Bunch current density, transverse
Jbr=fread(file,Nr*Nz,'float64');              
Jbr=reshape(Jbr,[Nz Nr]);
Jbr=Jbr';

% - Plasma current density, transverse
Jbckr=fread(file,Nr*Nz,'float64');              
Jbckr=reshape(Jbckr,[Nz Nr]);
Jbckr=Jbckr';

% - Bunch current density, longitudinal
Jbz=fread(file,Nr*Nz,'float64');              
Jbz=reshape(Jbz,[Nz Nr]);
Jbz=Jbz';

% - Plasma current density, longitudinal
Jbckz=fread(file,Nr*Nz,'float64');              
Jbckz=reshape(Jbckz,[Nz Nr]);
Jbckz=Jbckz';

% - Zstar
Zstar=fread(file,Nr*Nz,'float64');              
Zstar=reshape(Zstar,[Nz Nr]);
Zstar=Zstar';

% - rho_i
rho_i=fread(file,Nr*Nz,'float64');              
rho_i=reshape(rho_i,[Nz Nr]);
rho_i=rho_i';

fclose('all');
