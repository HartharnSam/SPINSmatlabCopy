function dd = cvdd(f,dim);
% Chebyshev derivative along dimension dim
if (~exist('dim')) dim=1; end;
sz = size(f);
perm = 1:length(sz);
perm(dim)=1;perm(1)=dim;
psz = sz(perm);

f = permute(f,perm);
f = reshape(f,[psz(1) prod(psz(2:length(psz)))]);

cc = cvt(real(f));

dd = zeros(size(cc));
N = psz(1);
dd(N-1,:) = (N-1)*cc(N,:);
for k = N-2:-1:1;
    dd(k,:) = dd(k+2,:) + 2*k*cc(k+1,:);
end

dd = cvt(real(dd))/(2*(N-1));
dd = reshape(dd,psz);
dd = permute(dd,perm);