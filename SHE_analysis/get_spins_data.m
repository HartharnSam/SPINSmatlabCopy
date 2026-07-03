function data = get_spins_data(var, ii, xminInd, xmaxInd, zInds)
%GET_SPINS_DATA - returns spins data, basically a wrapper for
%spins_reader_new, but it can do some computations for certain derived variables
if nargin < 3 || isempty(xminInd)
    xminInd = 1;
end
if nargin < 4 || isempty(xmaxInd)
    params = spins_params;
    xmaxInd = params.Nx;
end
if nargin < 5 || isempty(zInds)
    params = spins_params;
    zInds = 1:params.Nz;
end

switch lower(var)
    case 's'
        data = spins_reader_new('s', ii, xminInd:xmaxInd, zInds);
        data = data.*(data > 0);
    case 'enstrophy'
        try
            data = spins_reader_new('enst', ii, xminInd:xmaxInd, zInds);
        catch
            try
                data = 0.5*spins_reader_new('vorty', ii, xminInd:xmaxInd, zInds).^2;
            catch 
                spins_derivs('vorty', ii, true);
                data = 0.5*spins_reader_new('vorty', ii, xminInd:xmaxInd, zInds).^2;

            end
        end
    case 'ke'
        u = spins_reader_new('u', ii, xminInd:xmaxInd, zInds);
        w = spins_reader_new('w', ii, xminInd:xmaxInd, zInds);
        data = 0.5*(u.^2 + w.^2);
    case 'vorty'
        try
            data = spins_reader_new("vorty", ii, xminInd:xmaxInd, zInds);
        catch
            data = spins_derivs('vorty', ii);
            warning("Using spins_derivs, no save")
            data = data(xminInd:xmaxInd, zInds);
        end
    case 'speed'
        u = spins_reader_new('u',ii, xminInd:xmaxInd, zInds);
        w = spins_reader_new('w',ii, xminInd:xmaxInd, zInds);
        data = sqrt(u.^2 + w.^2);
    case 'diss'
        try
            data = spins_reader_new('diss', ii, xminInd:xmaxInd, zInds);
        catch 
            data = spins_derivs('diss', ii, true);
            data = data(xminInd:xmaxInd, zInds);
        end
        data = log10(data);

    case 'rho'
        try
            data = spins_reader_new('rho', ii, xminInd:xmaxInd, zInds);
        catch
            %rho0 = params.rho_0;
            data = (eqn_of_state(spins_reader_new('t', ii, xminInd:xmaxInd,...
                zInds), 0));
        end
        %data1 = data1 + -0.5*params.delta_rho * tanh((z-params.rho_loc)/params.dz_rho);
    case 'rho_z2'
        try
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, zInds).^2;
        catch
            spins_derivs('rho_z', ii, true);
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, zInds).^2;
        end
    case 'rho_zz2'
        try
            data = log10(spins_reader_new('rho_zz', ii, xminInd:xmaxInd, zInds).^2);
        catch
            spins_derivs('rho_zz', ii, true);
            data = log10(spins_reader_new('rho_zz', ii, xminInd:xmaxInd, zInds).^2);
        end
    case 'u_z2'
        try
            data = spins_reader_new('u_z', ii, xminInd:xmaxInd, zInds).^2;
        catch
            spins_derivs('u_z', ii, true);
            data = spins_reader_new('u_z', ii, xminInd:xmaxInd, zInds).^2;
            delete("u_z."+ii);
        end
    case 'grad_rho'
        try
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, zInds).^2;
            data = data + spins_reader_new('rho_x', ii, xminInd:xmaxInd, zInds).^2;
            data = sqrt(data);
        catch
            spins_derivs('rho_z', ii, true); spins_derivs('rho_x', ii, true);
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, zInds).^2;
            data = data + spins_reader_new('rho_x', ii, xminInd:xmaxInd, zInds).^2;
            data = sqrt(data);
            delete("rho_z."+ii);
        end
    otherwise
        try
            data = spins_reader_new(var, ii, xminInd:xmaxInd, zInds);
        catch
            try
                data = spins_derivs(var, ii, false, xminInd:xmaxInd, zInds);
            catch
                error(var+' not configured');
            end
        end
end
end