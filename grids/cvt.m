function cc = cvt(f)
% DCT-00 transform of input signal f, for Chebyshev transform
sz = size(f);
N=sz(1);
cc = fft([f; flipdim(f(2:N-1,:,:,:,:,:),1)]);
cc(N+(1:(N-2)),:,:,:,:,:)=[];
if isreal(f), cc = real(cc); end