% Input [1] 'Distance'
% Input [1] 'Data path'
% Reads the Phase space and 2D outputs of Architect
% Recommended use: Read_Architect_bin(Distance,pwd)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : F. Massimo, A. Marocchino
% Purpose       : read binary output from Architect
% Last modified : 17/3/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 function Read_Architect_bin(Distance,Folder)
 
 
filename_2D   = [Folder '\' num2str(Distance,'%0.6d') '_um.bin']; 
filename_PS   = [Folder '\' 'b_' num2str(Distance,'%0.6d') '_um.bin'];
filename_PS_r = [Folder '\' 'PS_' num2str(Distance,'%0.6d') '_um.bin'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Open and read 2D file              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = fopen(filename_2D,'rb','l');

% - domain axes
Nx=fread(file,1,'int32');                       % half-plane axis dim
Nz=fread(file,1,'int32');                       % z dim
x_axis_positive =fread(file,Nx/2,'float32'); 
z_axis          =fread(file,Nz,'float32');

% - bunch density
r=fread(file,Nx*Nz,'float32');              
r=reshape(r,[Nz Nx]);
r=r';

% - plasma electron density
n=fread(file,Nx*Nz,'float32');              
n=reshape(n,[Nz Nx]);
n=n';

% - Electric field, transverse
Er=fread(file,Nx*Nz,'float32');              
Er=reshape(Er,[Nz Nx]);
Er=Er';

% - Electric field, longitudinal
Ez=fread(file,Nx*Nz,'float32');              
Ez=reshape(Ez,[Nz Nx]);
Ez=Ez';

% - Magnetic field, azimuthal
Bphi=fread(file,Nx*Nz,'float32');              
Bphi=reshape(Bphi,[Nz Nx]);
Bphi=Bphi';

% - Bunch current density, transverse
Jbr=fread(file,Nx*Nz,'float32');              
Jbr=reshape(Jbr,[Nz Nx]);
Jbr=Jbr';

% - Plasma current density, transverse
Jer=fread(file,Nx*Nz,'float32');              
Jer=reshape(Jer,[Nz Nx]);
Jer=Jer';

% - Bunch current density, longitudinal
Jbz=fread(file,Nx*Nz,'float32');              
Jbz=reshape(Jbz,[Nz Nx]);
Jbz=Jbz';

% - Plasma current density, longitudinal
Jez=fread(file,Nx*Nz,'float32');              
Jez=reshape(Jez,[Nz Nx]);
Jez=Jez';

fclose(file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Open and read Phase Space files         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------
% ----------- Whole Phase Space --------------
% --------------------------------------------
file = fopen(filename_PS,'rb','l');

% --- Number of particles 
Np=fread(file,1,'int32');                       % whole phase space

% --- Reads phase space
Ps=fread(file,Np*9,'float32');              
Ps=reshape(Ps,[9 Np]);
Ps=Ps';

fclose(file);


% --------------------------------------------
% ---------- Reduced Phase Space -------------
% --------------------------------------------
file = fopen(filename_PS_r,'rb','l');
% --- Number of particles 
Npr=fread(file,1,'int32');                       % reduced phase space

% --- Reads reduced phase space
Ps_r=fread(file,Np*7,'float32');              
Ps_r=reshape(Ps_r,[7 Np]);
Ps_r=Ps_r';

fclose(file);

 end

