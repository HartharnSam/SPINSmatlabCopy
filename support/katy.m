function c = katy(N)
% Katy Perry colormap -- it's Hot and it's cold 
% (less middle-range contrast than darkjet)
%
% Courtesy John Ladan, 2015

if nargin < 1
   N = size(get(gcf,'colormap'),1);
end
h = hot(N);
c = flipud(fliplr(h));
h = padarray(h,[round(N/16*15),0],'pre');
c = padarray(c,[round(N/16*15),0],'post');
c = c+h;
c = c(1+round(N/8):end-round(N/8),:);
