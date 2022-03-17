function dd = fft_dx(f)
% Takes the derivative along the z (second) dimension
% via the Fourier transform, of the two-dimensional
% array f.  This function does -not- normalize by
% grid length, so a grid of x-length L should be
% multiplied by (2*pi)/L

[Nx Nz] = size(f);
persistent ks
% Create an array for the wavenumbers ks.  This is
% declared 'persistent' so it is -not- regenerated
% at each call.  In MATLAB, 'persistent' variables
% stay around after execution like 'global' ones, but
% do not show up in the global namespace.
if (~exist('ks') || size(ks,1) ~= Nx || size(ks,2) ~= Nz)
    ks = [0:floor((Nx-1)/2) -floor(Nx/2):-1]';
    ks = repmat(ks,[1 Nz]);
end

% Now, differentiation is a simple formula:

dd = real(ifft(1i.*ks.*fft(f,[],1),[],1));