function cd = cvd(f)
% Chebyshev derivative of function f

cc = cvt(real(f)); % Compute chebyshev / dct00 transform
cc = cc .* (abs(cc) > 1e-14*max(abs(cc(:))));
cd = zeros(size(cc));
N = length(cc);
cd(N-1) = (N-1)*cc(N);
for k = N-2:-1:1;
    cd(k) = (cd(k+2) + 2*(k)*cc(k+1));%/(1+(k==1));
end
cd = cvt(cd)/(2*(N-1)); % Inverse transform
