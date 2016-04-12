% *** read Architect binary SECTION output ***
% 
% input[1]  -> name with full path
% output[1] -> the entire output
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : A. Marocchino
% Purpose       : read SECTION binary output from Architect: option 2 with
% splitted fields
% Last modified : 17/3/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist,Nr,Nz,r_mesh,z_mesh,rho_b,n_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,Jbr,Jbckr,Jbz,Jbckz] = architect_read_bin_section_v2(full_name)

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
Er_bck=fread(file,Nr*Nz,'float32');
Er_bck=reshape(Er_bck,[Nz Nr]);
Er_bck=Er_bck';
Er_b=fread(file,Nr*Nz,'float32');
Er_b=reshape(Er_b,[Nz Nr]);
Er_b=Er_b';

% - Electric field, longitudinal
Ez=fread(file,Nr*Nz,'float32');              
Ez=reshape(Ez,[Nz Nr]);
Ez=Ez';
Ez_bck=fread(file,Nr*Nz,'float32');              
Ez_bck=reshape(Ez_bck,[Nz Nr]);
Ez_bck=Ez_bck';
Ez_b=fread(file,Nr*Nz,'float32');              
Ez_b=reshape(Ez_b,[Nz Nr]);
Ez_b=Ez_b';

% - Magnetic field, azimuthal
Bphi=fread(file,Nr*Nz,'float32');              
Bphi=reshape(Bphi,[Nz Nr]);
Bphi=Bphi';
Bphi_bck=fread(file,Nr*Nz,'float32');              
Bphi_bck=reshape(Bphi_bck,[Nz Nr]);
Bphi_bck=Bphi_bck';
Bphi_b=fread(file,Nr*Nz,'float32');              
Bphi_b=reshape(Bphi_b,[Nz Nr]);
Bphi_b=Bphi_b';

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
