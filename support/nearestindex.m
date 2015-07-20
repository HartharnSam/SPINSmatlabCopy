function index = nearestindex(list,value)
% NEARESTINDEX   Find the index of list that is the closest to value.
%    List is assumed to be a monotonic.
%
%    ind = nearestindex(x1d, x_i)
%
%    David Deepwell, 2015.

[~,index] = min(abs(list-value));
