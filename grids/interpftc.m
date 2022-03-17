function y = interpftc(x,ny,dim)
%INTERPFTC 1-D cell-centered interpolation using the FFT method
%   This function behaves as INTERFPT, save that it also applies
%   a phase shift such that the output grid is cell-centered -- i.e.
%   the leftmost boundary is one half-cell to the left of x(1)/y(1), and
%   the rightmoust boundary is one half-cell to the right of x(end)/y(end).
%   This solves an issue with interft such that an even/odd function of
%   size was not necessarily having the symmetry preserved upon 
%   interpolation

% From interpft -- permute the array as appropriate for multidimensional
% inputs, so that for this code we are only operating on the first
% dimension
if nargin==2,
  dim = find(size(x)>1,1,'first');
  if (isempty(dim)) dim = 1; end;
  %[x,nshifts] = shiftdim(x);
  %if isscalar(x), nshifts = 1; end % Return a row for a scalar
end
real_x = isreal(x);
%elseif nargin==3,
  perm = [dim:max(length(size(x)),dim) 1:dim-1];
  x = permute(x,perm);
%end

siz = size(x); %% Original size
[m,n] = size(x);  %% Size in a 2D-0

fx = fft(x,[],1); % Take Fourier transform of x
%kx = [0:floor(m/2) -floor((m-1)/2):-1]*2*pi; kx = kx(:);
% Shift the spectrum of fx, to set the first point as the left boundary
phs_fwd = -1/m/2; % Amount of phase shift needed on the forward transform

% If x has a nyquist frequency, divide it by two
if (~mod(m,2))
    fx(m/2+1,:) = fx(m/2+1,:)/2;
end

fy = zeros(ny,n);
if (ny <= m) % Downsampling
    fy = fx(1+[0:floor(ny/2) m+(-floor((ny-1)/2):-1)],:);
    if (~mod(ny,2)) % If there is a nyquist frequency on the grid
        % there is no unique choice for whether to use the positive
        % or negative-frequency equivalent, so choose the average.
        % Coincidentally, summing the two will also give the appropriate
        % scaling
        fy(1+ny/2,:) = fx(1+ny/2,:) + fx(m-ny/2+1,:);
    end
else % Upsampling
    % Positive frequencies
    fy(1:ceil((m+1)/2),:) = fx(1:ceil((m+1)/2),:);
    % Negative frequencies
    fy((ny-floor((m-1)/2)):ny,:) = fx((m-floor((m-1)/2)):m,:);
end

% Now, apply the proper phase shift to adjust the grid
phs_bk = 1/ny/2;
ky = [0:floor(ny/2) -floor((ny-1)/2):-1]*2*pi; ky = ky(:);
fy = full(spdiags(exp(1i*ky*(phs_fwd+phs_bk)),0,ny,ny)*fy);

if (~mod(ny,2) && real_x)
    % Adjust the nyquist frequency after the phase shift to ensure
    % a real output
    fy(ny/2+1,:) = real(fy(ny/2+1,:));
end

y = ifft(fy,[],1)*ny/m;

y = reshape(y,[ny siz(2:end)]);
    


% From interpft -- undo permutation done by multidimensional input
%if nargin==2,
%  y = reshape(y,[ones(1,nshifts) size(y,1) siz(2:end)]);
%elseif nargin==3,
  y = ipermute(y,perm);
%end