 function make_deriv_job(Directory, timestep)
cd(Directory)
params = spins_params; 
copyfile(['*.' num2str(timestep)], '../../TEMP_Ris');
copyfile('*grid', '../../TEMP_Ris');

orig_dir = cd;
cd('../../TEMP_Ris');
pause(1);
%% Read in SPINS.CONF_DERIV file line by line
fid = fopen('spins.conf_deriv', 'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

%% Set spatial parameters
A{4} = ['Lx = ' num2str(params.Lx)];
A{5} = ['Ly = ' num2str(params.Ly)];
A{6} = ['Lz = ' num2str(params.Lz)];
A{7} = ['Nx = ' num2str(params.Nx)];
A{8} = ['Ny = ' num2str(params.Ny)];
A{9} = ['Nz = ' num2str(params.Nz)];

%% Set some other things
A{20} = ['rho_0 = ' num2str(params.rho_0)];
A{24} = ['start_sequence = ' num2str(timestep)];
A{25} = ['final_sequence = ' num2str(timestep)];

% Set which derivatives to use
A{23} = 'deriv_files = u rho';
A{27} = 'deriv_x = false';
A{28} = 'deriv_y = false';
A{29} = 'deriv_z = true';
A{30} = 'do_vor_x = false';
A{31} = 'do_vor_y = false';
A{32} = 'do_vor_z = false';
A{33} = 'do_enstrophy = false';
A{34} = 'do_dissipation = true';
A{35} = 'do_vort_stretch = false';
A{36} = 'do_enst_stretch = false';

fid = fopen('spins.conf_deriv', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid, '%s', A{i});
        break
    else
        fprintf(fid, '%s \n', A{i});
    end
end

fclose(fid);
%%
disp('DO ROCKET THINGS')
pause
%%
movefile(['*.' num2str(timestep)], orig_dir)

cd(orig_dir);