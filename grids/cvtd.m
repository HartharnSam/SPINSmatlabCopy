function cc = cvtd(f,dim)
% DCT-00 transform of input signal f, for Chebyshev transform
if (~ exist('dim'))
    dim = 1;
end
cc = zeros(size(f));
sz = size(f);
numdim = length(sz);

perm = 1:numdim; % Compute permutation array
perm(1) = dim;
perm(dim) = 1; % Swap dimension 1 and dim

cc = permute(cvt(permute(f,perm)),perm);