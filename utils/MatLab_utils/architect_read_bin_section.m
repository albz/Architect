% *** read Architect binary SECTION output ***
% 
% input[1]  -> name with full path
% output[1] -> the entire output
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : A. Marocchino, F. Massimo
% Purpose       : read SECTION binary output from Architect
% Last modified : 17/3/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist,Nr,Nz,r_mesh,z_mesh,rho_b,n_bck,Er,Ez,Bphi,Jbr,Jbckr,Jbz,Jbckz] = architect_read_bin_section(full_name)

file = fopen(full_name,'rb','l');

% --- Output version
output_version=fread(file,1,'int32');

% --- Traveled distance
dist=fread(file,1,'int32');

% - domain axes
Nr=fread(file,1,'int32');                       % half-plane axis dim
Nz=fread(file,1,'int32');                       % z dim
r_mesh =fread(file,Nr,'float32');
z_mesh =fread(file,Nz,'float32');

% - bunch density
rho_b=fread(file,Nr*Nz,'float32');              
rho_b=reshape(rho_b,[Nz Nr]);
rho_b=rho_b';

% - plasma electron density
n_bck=fread(file,Nr*Nz,'float32');              
n_bck=reshape(n_bck,[Nz Nr]);
n_bck=n_bck';

% - Electric field, transverse
Er=fread(file,Nr*Nz,'float32');              
Er=reshape(Er,[Nz Nr]);
Er=Er';

% - Electric field, longitudinal
Ez=fread(file,Nr*Nz,'float32');              
Ez=reshape(Ez,[Nz Nr]);
Ez=Ez';

% - Magnetic field, azimuthal
Bphi=fread(file,Nr*Nz,'float32');              
Bphi=reshape(Bphi,[Nz Nr]);
Bphi=Bphi';

% - Bunch current density, transverse
Jbr=fread(file,Nr*Nz,'float32');              
Jbr=reshape(Jbr,[Nz Nr]);
Jbr=Jbr';

% - Plasma current density, transverse
Jbckr=fread(file,Nr*Nz,'float32');              
Jbckr=reshape(Jbckr,[Nz Nr]);
Jbckr=Jbckr';

% - Bunch current density, longitudinal
Jbz=fread(file,Nr*Nz,'float32');              
Jbz=reshape(Jbz,[Nz Nr]);
Jbz=Jbz';

% - Plasma current density, longitudinal
Jbckz=fread(file,Nr*Nz,'float32');              
Jbckz=reshape(Jbckz,[Nz Nr]);
Jbckz=Jbckz';

fclose('all');
