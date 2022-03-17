function dd = fft3_dy(f)
% Takes the derivative of three-dimensional array f via
% the FFT.  This does not normalize for domain length,
% so the result should be multiplied by 2*pi/Lx

[Nx Ny Nz] = size(f);
persistent ks
% Create an array for the wavenumbers ks.  This is
% declared 'persistent' so it is -not- regenerated
% at each call.  In MATLAB, 'persistent' variables
% stay around after execution like 'global' ones, but
% do not show up in the global namespace.
%if (~exist('ks') || size(ks,1) ~= Nx || size(ks,2) ~= Ny || size(ks,3) ~= Nz)
if (~exist('ks') || size(ks,2) ~= Ny)
    ks = [0:floor((Ny-1)/2) -floor(Ny/2):-1];
%    ks = repmat(ks,[Nx 1 Nz]);
end

% Now, differentiation is a simple formula:

%dd = real(ifft(1i.*ks.*fft(f,[],2),[],2));
dd = real(ifft(1i.*bsxfun(@times,ks,fft(f,[],2)),[],2));
