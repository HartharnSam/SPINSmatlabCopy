function [ape_tot] = calc_energetics(ii, xlimits, filter_rho, isCompare, ~)
% Script which employs the Winters (1995)-like sorting algorithm for a
% mapped case to calculate APE & BPE, and then also calculates KE and
% dissipation. Packages nicely into a structure
% Author: Marek Stastna

% Read in and cut down the grid
[x, z] = spinsgrid2d;
xind1 = nearest_index(x(:, 1), xlimits(1));
xind2 = nearest_index(x(:, 1), xlimits(2));
xinds = xind1:xind2;
x = x(xinds, :); z = z(xinds, :);

% Read in spins.conf - and calculate the related infos
params = spins_params;
sz = size(x);
Nx = sz(1);
Nz = sz(2);
Nzc = Nz - 1; % Nz for chebyshev
dx = x(2) - x(1);

% Read in rho, u and w data
%rho = rho_converter(spins_reader_new('rho',ii, xinds, []));
%rho = params.rho_0*(1+spins_reader_new('rho', ii, xinds, []));
rho = spins_reader_new('rho', ii, xinds, []);

u = spins_reader_new('u',ii, xinds, []);
w = spins_reader_new('w',ii, xinds, []);
rho_0 = params.rho_0;
% Filter out erroneously high density (that results from model
% mixing/filter)
if filter_rho
    myfact = 0.01; % percentage factor to allow overshoot by
    rho0 = (spins_reader_new('rho',0));
    rhomax = max(rho0(:));
    rhomin = min(rho0(:));
    delrho = rhomax-rhomin;
    %  Produces a "myind" - matrix of where density overshoots
    myind = 1.0*(rho>rhomax+(myfact*delrho))+1.0*(rho<rhomin-(myfact*delrho));
    if max(myind, [], 'all') == 1 % if there are overshoots
        myind = logical(myind);
        % Then remove them from ALL data
        u(myind) = NaN;
        w(myind) = NaN;
        rho(myind) = NaN;
        
    end
end

ke = 0.5*rho_0.*(u.^2+w.^2); % calculate Kinetic Energy

%% Calculate Chebyschev volumes
% Compute the area associated with each Chebyshev point using the values
% halfway between the point below and above
[~,z1dc] = cheb(Nzc);
[~,wci] = clencurt(Nzc);

% A normalised height applied to each Chebyshev point
arc(1) = 0.5*(z1dc(1)-z1dc(2)); % .5 comes from clencurt returning a total weight of 2
arc(Nzc+1) = arc(1);
for jj=2:(Nzc) % Central differencing
    arc(jj) = 0.5*(z1dc(jj-1)-z1dc(jj))+0.5*(z1dc(jj)-z1dc(jj+1)); % d_z1dc
end

% Then normalise this by the change in depth at each x position
arcphys = NaN(Nx, Nz);
for jj = 1:Nx
    Lznow = max(z(jj, :))-min(z(jj, :)); 
    arcphys(jj, :) = arc*Lznow; % Thickness of each grid point
end

if filter_rho && max(myind, [], 'all') == 1
    arcphys(myind) = NaN;
    arcphysv = arcphys(~isnan(arcphys(:)));
else
    arcphysv = arcphys(:);
end

%% Now do the sorting algorithm
% For chebyshev it takes a bit of work
if filter_rho
    [rhosortedc, rhosortedci] = sort(rho(~isnan(rho(:))), 'descend');
else
    [rhosortedc, rhosortedci] = sort(rho(:), 'descend');
end
% Now create the zsortedc
zsortedc = NaN(1, length(arcphysv));
zsortedc(1) = arcphysv(rhosortedci(1));
for jj = 2:(length(arcphysv))  % For each grid point
    zsortedc(jj) = zsortedc(jj-1)+arcphysv(rhosortedci(jj)); % point weights sorted by density
end
zsortedc = zsortedc-min(zsortedc);
zsortedc = zsortedc/max(zsortedc);
%TODO - Fix from here onwards

% This assumes the maximum depth is at the left end point
% Create a physical version of the sorted grid
zmin=min(z(1,:));
zmax=max(z(1,:));
zsortedc=zmin+(zmax-zmin)*zsortedc;
z_displace = zsortedc(find(rhosortedc > rhosortedc(end)+(.019*.04), 1, 'last'))
%mean_density = mean(rhosortedc);
%z_displace = zsortedc(nearest_index(rhosortedc, mean_density))
zsortedc = zsortedc-z_displace;
zmin = zmin-z_displace;
zmax = zmax-z_displace;
% now create a working grid so the repeated interpolations don't take
% forever - decrease the resolutions of zsortedc/rhosortedc
zsortedworking = linspace(zmin,zmax,201);
rhosortedworking = interp1(zsortedc,rhosortedc,zsortedworking,'linear');


% Finally create the background density profile at each x value
% And get the total APE as well as the tot APE at each x value
% I compute the KE info as well
ape_tot = 0;
ke_tot = 0;
bpe_tot = 0;
rhob = NaN(Nx, Nz);
ape_hor = NaN(1, Nx);
pe_tot = 0;
pe_hor = NaN(Nx, Nz);
ke_hor = NaN(1, Nx);
bpe_hor = NaN(1, Nx);
zi = z-z_displace;
z = z-params.min_z;

rho = rho_0*(1+rho);
for jj=1:Nx

    zminnow = min(z(jj,:));
    zmaxnow = max(z(jj,:));
    rhob(jj,:) = interp1(zsortedworking,rhosortedworking, zi(jj,:), 'spline'); % background (sorted) stratification
    rhob(jj, :) = rho_0*(1+rhob(jj, :));
    
    % get the local chain rule expression - Weight of each grid point
    winow = wci*(zmaxnow-zminnow)*0.5*dx;
    
    % Get APE
    % integrate vertically
    dummy = 9.81*sum(winow.*(rho(jj,:) - rhob(jj,:)).*zi(jj,:));
    % notice that if the surface is at z=0 z will have negative values
    ape_hor(jj) = dummy;
    ape_tot = ape_tot+dummy;
    
        % Get BPE
    % integrate vertically
    
    dummy = 9.81*sum(winow.*(rhob(jj,:)).*z(jj,:));
    % notice that if the surface is at z=0 z will have negative values
    bpe_hor(jj) = dummy;
    bpe_tot = bpe_tot+dummy;
    
    % Get PE
    dummy = 9.81*sum(winow.*rho(jj, :).*z(jj,:));
    pe_hor(jj) = dummy;
    pe_tot = pe_tot+dummy;
    
    % Get KE
    % integrate vertically
    dummy = sum(winow.*ke(jj,:));
    ke_hor(jj) = dummy;
    ke_tot = ke_tot+dummy;
end

if isCompare
    load('all_diagnos', 'all_diagnos')
    t_ind = nearest_index(all_diagnos.EnergyBudget.Time, ii);
    APE_tot = all_diagnos.EnergyBudget.APE_tot(t_ind);
    KE_tot = all_diagnos.EnergyBudget.KE_tot(t_ind);
    BPE_tot = all_diagnos.diagnos.BPE_tot(t_ind);
    PE_tot = all_diagnos.diagnos.PE_tot(t_ind);
    
    fprintf('Online APE:         %6.4f \n', APE_tot);
    fprintf('Offline APE:         %6.4f \n', ape_tot);
    fprintf('APE error:         %6.1f %% \n', 100*(APE_tot-ape_tot)/APE_tot);
    fprintf('Test APE:          %6.4f \n', pe_tot-bpe_tot);
    fprintf('Test APE error:         %6.1f %% \n', 100*(APE_tot-(pe_tot-bpe_tot))/APE_tot);
    
    disp('-----')
    fprintf('Online KE:         %6.4f \n', KE_tot);
    fprintf('Offline KE:         %6.4f \n', ke_tot);
    fprintf('KE error:         %6.1f %% \n', 100*(KE_tot-ke_tot)/KE_tot);
    
    disp('-----')
    fprintf('Online BPE:         %6.4f \n', BPE_tot);
    fprintf('Offline BPE:         %6.4f \n', bpe_tot);
    fprintf('BPE error:         %6.1f %% \n', 100*(BPE_tot-bpe_tot)/BPE_tot);
    
    disp('-----')
    fprintf('Online PE:         %6.4f \n', PE_tot);
    fprintf('Offline PE:        %6.4f \n', pe_tot);
    fprintf('PE error:          %6.1f %% \n', 100*(PE_tot-pe_tot)/PE_tot);
    
end
ape_tot = struct('APE_Total', ape_tot, 'PE_Total', pe_tot, 'KE_Total', ke_tot, ...
    'BPE_Total', bpe_tot);

end
