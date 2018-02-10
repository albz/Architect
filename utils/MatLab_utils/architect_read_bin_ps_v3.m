% *** read Architect binary SECTION output ***
% 
% input[1]  -> name with full path
% output[1] -> the entire output
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : A. Marocchino, F. Massimo
% Purpose       : read PHASE SPACE binary output from Architect
% Last modified : 17/3/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist,x,y,z,px,py,pz,bunch_id,cut,dcut,bunch_charges,macro_particle_charge,macro_particle_nume] = architect_read_bin_ps_v3(full_name)
    file = fopen(full_name,'rb','l');
    

    % --- Output version
    output_version=fread(file,1,'int32');
    
    % --- Bunch number and Charges
    n_bunches=fread(file,1,'int32');
    bunch_charges=fread(file,n_bunches,'float64');
    
    % --- Traveled distance
    dist=fread(file,1,'int32');
    
    % --- Number of particles
    Np=fread(file,1,'int32');

    % --- Reads phase space
    ps=fread(file,Np*11,'float64');
    ps=reshape(ps,[11 Np]);
    
    x=ps(1,:);
    y=ps(2,:);
    z=ps(3,:);
    px=ps(4,:);
    py=ps(5,:);
    pz=ps(6,:);
    bunch_id=ps(7,:);
    cut=ps(8,:);
    dcut=ps(9,:);
    macro_particle_charge=ps(10,:);
    macro_particle_nume=ps(11,:);

    fclose('all');
