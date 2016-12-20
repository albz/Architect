% *** read Architect binary output ***
%
% input[1] -> 'PS' or 'section'
%             if missing: selecting 2D outputs
%
% input[2] -> number (distance travelled in um)
%             if missing: last availabe output
% input[3] -> path
%             if missing: pwd
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : F. Massimo, A. Marocchino
% Purpose       : read binary output from Architect
% Last modified : 17/3/2015
%               : 2015/03/17 -> Complete rewritten with more flexible structure
%               : 2015/03/17 -> reading sub-structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function architect_read_bin(kind,distance,path)
%---***---%
% nargin=2
if( nargin==0 )
    path  = '.';
    files = dir(fullfile(path, 'out', '2D', '*.arch'));
    count = size(files,1);
    name  = files(count).name;
    full_name = fullfile(path, 'out', '2D', name);
    kind  ='section';
elseif( nargin==1 )
    path = '.';
    if strcmp(kind,'section')
        files = dir(fullfile(path, 'out', '2D', '*.arch'));
        count = size(files,1);
        name      = files(count).name;
        full_name = fullfile(path, 'out', '2D', name);
    elseif strcmp(kind,'PS')
        files = dir(fullfile(path, 'out', 'PS', '*.arch'));
        count = size(files,1);
        name      = files(count).name;
        full_name = fullfile(path, 'out', 'PS', name);
    end
elseif( nargin==2 )
    path = '.';
    if strcmp(kind,'section')
        files = dir(fullfile(path, 'out', '2D', '*.arch'));
        name      = files(distance).name;
        full_name = fullfile(path, 'out', '2D', name);
    elseif strcmp(kind,'PS')
        files = dir(fullfile(path, 'out', 'PS', '*.arch'));
        name      = files(distance).name;
        full_name = fullfile(path, 'out', 'PS', name);
    end
elseif( nargin==3 )
    if strcmp(kind,'section')
        files = dir(fullfile(path, 'out', '2D', '*.arch'));
        name      = files(distance).name;
        full_name = fullfile(path, 'out', '2D', name);
    elseif strcmp(kind,'PS')
        files = dir(fullfile(path, 'out', 'PS', '*.arch'));
        name      = files(distance).name;
        full_name = fullfile(path, 'out', 'PS', name);
    end
end

fprintf('file name: %s \n',full_name);

%---***---%
% if(strcmp(kind,'section') && nargin >=2 )
%     name      = [num2str(distance, '%0.6i'),'_um.bin'];
%     full_name = fullfile(path, 'out', '2D', name);
% end
% if(strcmp(kind,'ps')  && nargin >=2 )
%     name      = ['PS_',num2str(distance, '%0.6i'),'_um.bin'];
%     full_name = fullfile(path, 'out', 'PS', name);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              section                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(kind,'section'))
    full_name
    file = fopen(full_name,'rb','l');
    % --- Output version
    output_version=fread(file,1,'int32');
    fclose('all');

    if(output_version == 1)
        [dist,Nr,Nz,r_mesh,z_mesh,rho_b,rho_bck,Er,Ez,Bphi,Jbr,Jbckr,Jbz,Jbckz] = architect_read_bin_section(full_name);

        assignin('base', 'Nr', Nr);
        assignin('base', 'Nz', Nz);
        assignin('base', 'r_mesh', r_mesh);
        assignin('base', 'z_mesh', z_mesh);
        assignin('base', 'rho_b', rho_b);
        assignin('base', 'rho_bck', rho_bck);
        assignin('base', 'Er', Er);
        assignin('base', 'Ez', Ez);
        assignin('base', 'Bphi', Bphi);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbckr', Jbckr);
        assignin('base', 'Jbz', Jbz);
        assignin('base', 'Jbckz', Jbckz);
        assignin('base', 'dist', dist);


    elseif (output_version==2)

        [dist,Nr,Nz,r_mesh,z_mesh,rho_b,rho_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,Jbr,Jbckr,Jbz,Jbckz] = architect_read_bin_section_v2(full_name);

        assignin('base', 'Nr', Nr);
        assignin('base', 'Nz', Nz);
        assignin('base', 'r_mesh', r_mesh);
        assignin('base', 'z_mesh', z_mesh);
        assignin('base', 'rho', rho_b+rho_bck);
        assignin('base', 'rho_b', rho_b);
        assignin('base', 'rho_bck', rho_bck);
        assignin('base', 'Er', Er);
        assignin('base', 'Er_bck', Er_bck);
        assignin('base', 'Er_b', Er_b);
        assignin('base', 'Ez', Ez);
        assignin('base', 'Ez_bck', Ez_bck);
        assignin('base', 'Ez_b', Ez_b);
        assignin('base', 'Bphi', Bphi);
        assignin('base', 'Bphi_bck', Bphi_bck);
        assignin('base', 'Bphi_b', Bphi_b);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbckr', Jbckr);
        assignin('base', 'Jbz', Jbz);
        assignin('base', 'Jbckz', Jbckz);
        assignin('base', 'dist', dist);



    elseif (output_version==3)

        [dist,Nr,Nz,r_mesh,z_mesh,rho_b,rho_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz] = architect_read_bin_section_v3(full_name);

        assignin('base', 'Nr', Nr);
        assignin('base', 'Nz', Nz);
        assignin('base', 'r_mesh', r_mesh);
        assignin('base', 'z_mesh', z_mesh);
        assignin('base', 'rho', rho_b+rho_bck);
        assignin('base', 'rho_b', rho_b);
        assignin('base', 'rho_bck', rho_bck);
        assignin('base', 'Er', Er);
        assignin('base', 'Er_bck', Er_bck);
        assignin('base', 'Er_b', Er_b);
        assignin('base', 'Ez', Ez);
        assignin('base', 'Ez_bck', Ez_bck);
        assignin('base', 'Ez_b', Ez_b);
        assignin('base', 'Bphi', Bphi);
        assignin('base', 'Bphi_bck', Bphi_bck);
        assignin('base', 'Bphi_b', Bphi_b);
        assignin('base', 'B_ex_poloidal', B_ex_poloidal);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbckr', Jbckr);
        assignin('base', 'Jbz', Jbz);
        assignin('base', 'Jbckz', Jbckz);
        assignin('base', 'dist', dist);

    
    
    elseif (output_version==4)

        [dist,Nr,Nz,r_mesh,z_mesh,rho_b,rho_bck,Er,Er_bck,Er_b,Ez,Ez_bck,Ez_b,Bphi,Bphi_bck,Bphi_b,B_ex_poloidal,Jbr,Jbckr,Jbz,Jbckz,Zstar,rho_i] = architect_read_bin_section_v4(full_name);

        assignin('base', 'Nr', Nr);
        assignin('base', 'Nz', Nz);
        assignin('base', 'r_mesh', r_mesh);
        assignin('base', 'z_mesh', z_mesh);
        assignin('base', 'rho', rho_b+rho_bck);
        assignin('base', 'rho_b', rho_b);
        assignin('base', 'rho_bck', rho_bck);
        assignin('base', 'Er', Er);
        assignin('base', 'Er_bck', Er_bck);
        assignin('base', 'Er_b', Er_b);
        assignin('base', 'Ez', Ez);
        assignin('base', 'Ez_bck', Ez_bck);
        assignin('base', 'Ez_b', Ez_b);
        assignin('base', 'Bphi', Bphi);
        assignin('base', 'Bphi_bck', Bphi_bck);
        assignin('base', 'Bphi_b', Bphi_b);
        assignin('base', 'B_ex_poloidal', B_ex_poloidal);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbr', Jbr);
        assignin('base', 'Jbckr', Jbckr);
        assignin('base', 'Jbz', Jbz);
        assignin('base', 'Jbckz', Jbckz);
        assignin('base', 'Zstar', Zstar);
        assignin('base', 'rho_i', rho_i);
        assignin('base', 'dist', dist);        
    
    end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Phase Space                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(kind,'PS'))
    full_name
    file = fopen(full_name,'rb','l');
    % --- Output version
    output_version=fread(file,1,'int32');
    fclose('all');
    if(output_version == 1) [dist,x,y,z,px,py,pz,bunch_id,cut,dcut] = architect_read_bin_ps(full_name); end
    if(output_version == 2) [dist,x,y,z,px,py,pz,bunch_id,cut,dcut,bunch_charges] = architect_read_bin_ps_v2(full_name); end

    assignin('base', 'dist', dist);
    assignin('base', 'x', x);
    assignin('base', 'y', y);
    assignin('base', 'z', z);
    assignin('base', 'px', px);
    assignin('base', 'py', py);
    assignin('base', 'pz', pz);
    assignin('base', 'bunch_id', bunch_id);
    assignin('base', 'cut', cut);
    assignin('base', 'dcut', dcut);
    if(output_version == 2) assignin('base', 'bunch_charges', bunch_charges); end
end
