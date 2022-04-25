function [KdV, c_0] = calc_kdv_depthchange(depths)
%CALC_KDV_DEPTHCHANGE - % pseudospectral solution of the vertical structure
% of linear internal waves:
%  lambda(-D^2+k^2 I)phi = N^2(z)phi
% on 0<z<H with phi(0)=phi(H)=0
% This script varies the total depth for a fixed stratification
%
% Inputs:
%    depths - Vector of depths to calculate KdV for
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Marek Stastna
% Department of Applied Mathematics, University of Waterloo
% email address: mmstastna@uwaterloo.ca
% GitHub: https://github.com/HartharnSam
% 22-Mar-2022; Last revision: 22-Mar-2022
% MATLAB Version: 9.12.0.1884302 (R2022a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clf
format long, format compact

% set this to 1 if you want live plots
PLOT_NOW=0;
% This is what you need to get the differentiation matrix and the grid
N=200;
params = spins_params;
Hdeep = params.Lz;
if nargin < 1
    depths = Hdeep;
end
g = params.g;
delrho = params.delta_rho/2;
% z0deep = 0.3; % z0 - HAB of pycnocline at full depth
z0deep = -(params.min_z - params.pyc_loc);
%d = 0.01; % d - pycnocline halfthickness
d = params.h_halfwidth;

[D,zc] = cheb(N);
D2 = D^2; D2 = D2(2:N,2:N);

[~, w]=clencurt(N);
% here are the physical parameters and the scaling to the computational
% domain [-1, 1]
Hs = depths;
mylen = length(Hs);

%Hs = linspace(Hdeep,0.2*Hdeep,mylen);
%% Pre-allocate
r10_1 = NaN(1, mylen); r01_1 = NaN(1, mylen);
r10_2 = NaN(1, mylen); r01_2 = NaN(1, mylen);
c1s = r10_2; c2s = c1s;

%% Run for each depth
for cntr=1:mylen
    H=Hs(cntr);
    dzcdzp=2/H;
    dzpdzc=(H/2);
    z0=z0deep+(H-Hdeep);
    zphys=0.5*H*(zc+1);
    %dz_num = 1e-8*H;

    % inline functions for the density and the derivative of the density
    myrho=@(z) 1-delrho*tanh((z-z0)/d);
    myn2=@(z) (g*delrho/d)*sech((z-z0)/d).^2;
    if cntr == 1
        figure
        plot(myrho(zphys), zphys);
    end
    n2physical=myn2(zphys);
    %n2max=max(n2physical);

    % make up the matrices for the e-val prog.

    %define B
    B = -D2*(2/H)^2;
    %define A
    A = diag(n2physical(2:end-1));
    % Solve the e-val prob
    [ev, ee] = eig(A,B);
    [cs, csi]=sort(sqrt(diag(ee)),'descend');
    c1=cs(1);c2=cs(2);%c3=cs(3);
    % This makes sure that the eigenfunction has a maximum of 1
    phi1=ev(:,csi(1));
    %mxphi1=max(phi1);
    mnphi1=min(phi1);
    mxabs=max(abs(phi1));
    if abs(mnphi1)==mxabs
        phi1=-phi1/mxabs;
    else
        phi1=phi1/mxabs;
    end
    phi2=ev(:,csi(2));
    %mxphi2=max(phi2);
    mnphi2=min(phi2);
    mxabs2=max(abs(phi2));
    if abs(mnphi2)==mxabs2
        phi2=-phi2/mxabs2;
    else
        phi2=phi2/mxabs2;
    end
    if PLOT_NOW==1
        figure
        betterplots
        subplot(1,2,1)
        plot(n2physical,zphys+(Hdeep-H)),grid on,ylabel('z (m)'),xlabel('N^2')
        subplot(1,2,2)
        plot(phi1,zphys(2:end-1)+(Hdeep-H)),grid on,ylabel('z (m)'),xlabel('\phi')
        title(['c1 = ' num2str(c1,4)])
    end
    % Here is WNL stuff
    phi1p = D*[0;phi1;0]*dzcdzp;
    S1 = sum(w'.*(phi1p.^2)*dzpdzc);
    r10_1(cntr) = -0.75*sum(w'.*(phi1p.^3)*dzpdzc)/S1;
    r01_1(cntr) = -0.5*c1*sum(w'.*([0;phi1;0].^2)*dzpdzc)/S1;
    phi2p = D*[0;phi2;0]*dzcdzp;
    S2 = sum(w'.*(phi2p.^2)*dzpdzc);
    r10_2(cntr) = -0.75*sum(w'.*(phi2p.^3)*dzpdzc)/S1;
    r01_2(cntr) = -0.5*c2*sum(w'.*([0;phi2;0].^2)*dzpdzc)/S2;
    c1s(cntr) = c1;
    c2s(cntr) = c2;
end
beta = r01_2;
alpha = r10_2.*c1;

if mylen > 1
    figure(2)
    %betterplots;
    tiledlayout(3, 1);
    A1 = nexttile;
    plot(A1, Hs,c1s,'o-'),grid on,ylabel('c1 (m/s)'), xlabel('H (m)'); xlim([.15 .3]); hold on
    B2 = nexttile;
    plot(B2, Hs, alpha, 'o-'); xlim([.15 .3]); hold on;
    C3 = nexttile;
    plot(C3, Hs, beta, 'o-'); xlim([.15 .3]); hold on;
end
KdV = struct('Depths', Hs, 'c1s', c1s, 'c2s', c2s, 'r10_2', r10_2, 'r10_1', r10_1, 'r01_1', r01_1, 'r01_2', r01_2, 'c_0', c1s(1), 'alpha', alpha, 'beta', beta);
save('KdV', 'KdV');
fprintf('c_lw = %2.2f m/s \n', c1s(1));
