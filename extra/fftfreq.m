function ks = fftfreq(N, d, keep_kmax);
%  FFTFREQ  Return the Discrete Fourier Transform sample frequencies
%
%  Inputs:
%    'N'         - Number of points
%    'd'         - Sample spacing
%    'keep_kmax' - (Optional, boolean) Set extreme wavenumbers to zero
%
%  Outputs:
%    'ks'        - Vector of wavenumbers

% manage optional argument (keep_kmax)
if nargin < 3
    keep_kmax = true;
end

% minimum wavenumber
dk = 2*pi/(N*d);

if ~mod(N,2) % Even
    ks = [0:N/2 -N/2+1:-1]*dk;
else % Odd
    ks = [0:(N-1)/2 -(N-1)/2:-1]*dk;
end

if ~keep_kmax 
    % fix the extreme wavenumbers to be zero
    % to stop aliasing. See Trefethen
    if ~mod(N,2) % Even
        ks(N/2+1) = 0;
    else % Odd
        ks((N-1)/2+1:(N-1)/2+2) = [0 0];
    end
end
