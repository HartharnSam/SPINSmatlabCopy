function out = cvtinterp(f,N)
% Interpolates f to a grid of N points using chebyshev
% transforms

cout = zeros(N,1);
cout(1:length(f)) = cvt(f);
cout(length(f)) = cout(length(f))/2;
cout = cout(1:N);
cout(end) = cout(end)*2;
out = cvt(cout)/(2*(length(f)-1));

% Normalize here