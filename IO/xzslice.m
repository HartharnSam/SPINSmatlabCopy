function myslice = xzslice(f,jj)
% Returns the jjth xz slice of a 3D spins field
myslice=squeeze(f(:,jj,:));
