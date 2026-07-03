function data = spins_derivs(derivative, ii, save, xinds, zinds)
% Calculates derivatives offline. Currently only set up for 2D
% calculate dx dy dz as appropriate
% switch for derivatives (dissipation, vorticity, 
%
% Other m-files required: spins_params; spins_reader_new; get_grad2,
% spins_writer
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 12-Oct-2021; Last revision: 05-Sept-2024
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5
arguments
    derivative (1, 1) string 
    ii (1, 1) uint16 
    save (1, 1) logical = false
    xinds (:, 1) double = []
    zinds (:, 1) double = []
end
if isempty(save)
    save = false;
end
params = spins_params; 

switch lower(derivative)
    case 'ri' % Richardson Number
        g_rho0 = -params.g/params.rho_0;
        u = spins_reader_new('u', ii);    
        try
            rho = (spins_reader_new('rho', ii)+1)*params.rho_0;
        catch
            rho = nleos(spins_reader_new('s', ii), spins_reader_new('t', ii));
        end
        [~, du_dz] = get_grad2(u, xinds, zinds);
        [~, drho_dz] = get_grad2(rho, xinds, zinds);
        N_sq = g_rho0 * drho_dz;
        N_sq(abs(N_sq)<1e-2) = NaN; 
        du_dz(abs(du_dz)<21e-4) = NaN; 

        data = N_sq./(du_dz.^2);
        %data = smoothdata(data, 2, 'rlowess', 5);
    case 'n2'
        g_rho0 = -params.g/params.rho_0;
        rho = (spins_reader_new('rho', ii)+1)*params.rho_0;
        [~, drho_dz] = get_grad2(rho, xinds, zinds);
        data = g_rho0 * drho_dz;
        %data = smoothdata(data, 2, 'rlowess', 5);

    case 'vorty' % Vorticity
        u = spins_reader_new('u', ii);        
        w = spins_reader_new('w', ii);
        [~, du_dz] = get_grad2(u, xinds, zinds);
        [dw_dx] = get_grad2(w, xinds, zinds);
        data = du_dz-dw_dx;
        
    case 'diss' % Dissipation
        nu = params.visco*params.rho_0;
        u = spins_reader_new('u', ii);
        w = spins_reader_new('w', ii);
        [du_dx, du_dz] = get_grad2(u, xinds, zinds);
        [dw_dx, dw_dz] = get_grad2(w, xinds, zinds);
        data = 2*nu*(du_dx.^2 + dw_dz.^2 + 2*(.5*(du_dz + dw_dx)).^2);
    otherwise
        split_string = split(lower(derivative), '_');
        if ~(length(split_string)==2)
            error('derivative type not configured, try ri, vorty, diss, or var_x type')
        end

        raw_data = spins_reader_new(split_string{1}, ii);
        
        [~, z] = spinsgrid2d;
%        raw_data = raw_data + -0.5*params.delta_rho * tanh((z-params.rho_loc)/params.dz_rho);
        switch split_string{2}
            case 'x'
                data = get_grad2(raw_data, xinds, zinds);
            case 'z'
                [~, data] = get_grad2(raw_data, xinds, zinds);
            case 'zz'
                [~, data] = get_grad2(raw_data);
                [~, data] = get_grad2(data, xinds, zinds);
            case 'absgrad'
                [d_dx, d_dz] = get_grad2(raw_data, xinds, zinds);
                data = sqrt(d_dx.^2 + d_dz.^2);
            otherwise
                error('derivative type not configured, try ri, vorty, diss, or var_x type')
        end
end
% Save as SPINS output type file
filename = derivative+"."+ii;
if save
    spins_writer(filename, data);
end