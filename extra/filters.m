% define functions for different spectral filters

%% Hard cut-off filters
% k:  wavenumber
% kc: critical cut-off
f_step = @(kc,k) 1.*(abs(k)<=kc) + 0.*(abs(k)>kc);  % general cut-off
f_23   = @(k) f_step(2/3, k);                       % two-thirds cut-off

%% Gaussian filter
% a: filter strength
% b: filter order
f_gauss = @(a,b,kc,k) 1.*(abs(k)<=kc) + min(1, exp(-a*(abs(k)-kc).^b./(1-kc).^b)).*(abs(k)>kc);

%% hyperviscosity (a is adaptive in spins)
f_hyper = @(a,b,k) exp(-a*abs(k).^b);

%% Bump filters
f_bump  = @(a,b,kc,k) 1.*(abs(k)<=kc) + min(1, exp(a*(1-1./(1-((abs(k)-kc)./(1-kc)).^b)))).*(abs(k)>kc);
% k1: cut-off to begin filtering
% k2: cut-off above which nothing survives
f_bump2 = @(a,b,k1,k2,k) 1.*(abs(k)<=k1) ...
          + min(1, exp(a*(1-1./(1-((abs(k)-k1)./(k2-k1)).^b)))).*(abs(k)>k1 & abs(k)<k2) ...
          + 0.*(abs(k)>=k2);
