%Test script to demonstrate the different EOS functions available in
%SPINSmatlab
% Requires  nleos.m, quadeos.m, lineos.m, and compute_eos_coeffs.m
%
% AUTHOR:  Andrew Grace 2021-06-08  (andrew.grace@uwaterloo.ca)

close all;
%%
Tmin = 0;
Tmax = 40;
Smin = 0;
Smax = 45;


%% Fresh Comparison

                %Compares equations of state in fresh water. 
                %Plots nonlinear equation of state in fresh water over 0 to 40 deg C,
                %quadratic eos over 0 to 10 deg C, and the linear eos at two reference
                %temperatures. Reference salinity is set to 0. This is important because
                %if the reference salinity is anything else, the computed reference density will be
                %off.

figure(1);

T = linspace(Tmin,Tmax,4000);
S = 0;
T01 = 15;   %Reference temp 1
T02 = 30;   %Reference temp 2
S0 = 0;     %Reference Salinity

rho_lin1 = lineos(S,T,S0,T01);  %Compute linear eos at reference temp 1
rho_lin2 = lineos(S,T,S0,T02);  %Compute linear eos at reference temp 2
rho_quad = quadeos(T(1:1000));  %Compute quadratic eos over 0 to 10 deg C
rho_nl = nleos(S,T);            %Compute nonlinear eos for 0 salinity and temperature range
                                %If either S or T is a scalar and the other isn't, nleos.m assumes the
                                %scalar represents a uniform value and makes the sizing correct.

    
%plot    
plot(T(1:1000),rho_quad,'linew',2)
hold on
plot(T,rho_lin1,'linew',1);
plot(T,rho_lin2,'linew',1);
plot(T,rho_nl,'--k','linew',1);

%plot markers at reference temperatures
scatter([T01 T02],[nleos(S0,T01) nleos(S0,T02)],70,'ko','filled')

legend('Quadratic Fresh','Linear Fresh 1','Linear Fresh 2','Full Nonlinear',...
    'interpreter','latex');

grid on

xlabel('$T$ ($^\circ$C)','interpreter','latex');
ylabel('$\rho$ (kg$/$m$^3$)','interpreter','latex');

%% Uniform Salty, but variations in Temperature

            %The block of code below shows the behaviour of the functions when a
            %uniform but non-zero salinity is assumed. Since the full eos is nonlinear,
            %the reference temperature, expansion/contraction coefficients will be
            %different.

figure(2)
T = linspace(Tmin,Tmax,400);

T01 = 15;   %Reference temp 1
T02 = 30;   %Reference temp 2
S01 = 5;    %Reference salinity 1
S02 = 35;   %Reference salinity 2


%Compute linear eos at reference temps and salinities

rho_lin1 = lineos(S01,T,S01,T01);
rho_lin2 = lineos(S02,T,S02,T01);

rho_lin3 = lineos(S01,T,S01,T02);
rho_lin4 = lineos(S02,T,S02,T02);

%Compute nonlinear eos at reference temps and salinities

rho_nl1 = nleos(S01,T);
rho_nl2 = nleos(S02,T);

hold on
%plot
plot(T,rho_lin1,'-r','linew',1);
plot(T,rho_lin2,'--r','linew',1);

plot(T,rho_lin3,'-b','linew',1);
plot(T,rho_lin4,'--b','linew',1);

plot(T,rho_nl1,'-k','linew',1);
plot(T,rho_nl2,'--k','linew',1);

%plot markers at reference temperatures
scatter([T01 T01],[nleos(S01,T01) nleos(S02,T01)],70,'ko','filled')
scatter([T02 T02],[nleos(S01,T02) nleos(S02,T02)],70,'ko','filled')

legend('Linear 1','Linear 2','Linear 3','Linear 4',...
    'Nonlinear 1','Nonlinear 2','interpreter','latex');

grid on
xlabel('$T$ ($^\circ$C)','interpreter','latex');
ylabel('$\rho$ (kg$/$m$^3$)','interpreter','latex');

%% Salty 

            %The block of code below shows the behaviour for variations in salinity and
            %otherwise uniform temperature.

figure(3) 

T01 = 8;    %Reference temp 1
T02 = 25;   %Reference temp 2
S = linspace(Smin,Smax,450);

S01 = 5;    %Reference Salinity 1
S02 = 35;   %Reference Salinity 2

%Compute linear eos at reference temps and salinities

rho_lin1 = lineos(S,T01,S01,T01);
rho_lin2 = lineos(S,T01,S02,T01);
rho_lin3 = lineos(S,T02,S01,T02);
rho_lin4 = lineos(S,T02,S02,T02);


%Compute nonlinear eos at reference temps

rho_nl1 = nleos(S,T01);
rho_nl2 = nleos(S,T02);

%plot
plot(S,rho_lin1,'--r','linew',1); hold on
plot(S,rho_lin2,'--b','linew',1);
plot(S,rho_lin3,'-r','linew',1);
plot(S,rho_lin4,'-b','linew',1);

plot(S,rho_nl1,'--k','linew',1); hold on
plot(S,rho_nl2,'-k','linew',1);

%plot markers at reference salinities
scatter([S01 S02],[nleos(S01,T01) nleos(S02,T01)],70,'ko','filled')
scatter([S01 S02],[nleos(S01,T02) nleos(S02,T02)],70,'ko','filled')

legend('Linear 1','Linear 2','Linear 3','Linear 4',...
    'Nonlinear 1','Nonlinear 2','interpreter','latex');

grid on
xlabel('$S$ (psu)','interpreter','latex');
ylabel('$\rho$ (kg$/$m$^3$)','interpreter','latex');

%% Full NLEOS 

                %Plots a contour plot with labels for full nonlinear eos

[ss,tt] = meshgrid(S,T);

figure(4)

[C,h] = contour(tt,ss,nleos(ss,tt),20);
clabel(C,h,'LabelSpacing',800);

xlabel('$T$ ($^\circ$C)','interpreter','latex');
ylabel('$S$ (psu)','interpreter','latex');
