function [dx_out dy_out dz_out] = get_grad(in_func)
% Compute the gradient (f_x, f_y, f_z) of in_func, using appropriate
% mappings, on a fourier/fourier (unmapped)/cheby grid

persistent Jax Jbx Jaz Jbz Lx Ly Lz Px Py Pz
if (isempty(Jax) || any(size(Jax) ~= size(in_func)))

    % Compute Jacobians
    xgrid = xgrid_reader();
    if ndims(xgrid) >2
        [Nx Ny Nz] = size(xgrid);
        ygrid = ygrid_reader();

    else
        [Nx Nz] = size(xgrid);
                ygrid = [1 1];
        Ny = 1;
    end

    zgrid = zgrid_reader();
    
    
    % z includes boundaries, so Lz is the max - min
    Lz = max(zgrid(:))-min(zgrid(:));
    % Check whether the grid is increasing or decreasing
    Pz = sign(zgrid(1,1,1)-zgrid(1,1,end));
    
    % x and y do not, so there are really N-1 intervals;
    % Lx and Ly need scaled appropriately
    Lx = (max(xgrid(:)) - min(xgrid(:)))*(Nx)/(Nx-1);
    Px = sign(xgrid(end,1,1)-xgrid(1,1,1));
    if Ny>1
        Ly = (max(ygrid(:)) - min(ygrid(:)))*(Ny)/(Ny-1);
    else
        Ly = 1;
    end
    Py = sign(ygrid(1,end,1)-ygrid(1,1,1));
    
    % Now, create a temporary unmapped grid in x; this is necessary
    % because a linear slope (0 1 2 3) isn't differentiable via FFT
    % because of the implicit discontinuity at the boundary.  Instead
    % since we know the grid is actually increasing, we want to
    % explicitly remove that factor.
    
    Cgx = repmat(Px*Lx*(1:Nx)'/Nx,[1 Ny Nz]);
    
    % Now, compute dx/da -- derivative of x wrt the first coordinate
    dxda = fft3_dx(xgrid-Cgx) + Px*Lx/2/pi;
    % and dx/db
    dxdb = cvdd(xgrid,3);
    
    % Also, dz/da -- derivative of z wrt coordinate 1
    dzda = fft3_dx(zgrid);
    dzdb = cvdd(zgrid,3);
    
    % [dxda dxdb; dzda dzdb] forms a 2x2 matrix; the inverse of this
    % is the Jacobian transform of interest.  This is invertible, see
    % (2.74) of my thesis.
    
    JJ = dxda.*dzdb - dxdb.*dzda;
    Jax = dzdb./JJ;
    Jbx = -dzda./JJ;
    Jaz = -dxdb./JJ;
    Jbz = dxda./JJ;
end

% With that mess out of the way, compute derivatives

da = fft3_dx(in_func);
dy = fft3_dy(in_func);
db = cvdd(in_func,3);

dx_out = da.*Jax + db.*Jbx;
dz_out = da.*Jaz + db.*Jbz;
dy_out = dy.*(2*pi/Ly*Py);
    