function index = nearest_index(list, value)
%  NEARESTINDEX   Find the index of 'list' that is the closest to 'value'.
%
%  Usage:
%    ind = nearest_index(x1d, x_i)
%
%  Inputs:
%    'list'     - a monotonic vector (or any higher dimensional matrix)
%    'value'	- a number
%
%  Outputs:
%    'index'	- index of 'list' who's value is closest to 'value'
%
%  David Deepwell, 2015

    [~,index] = min(abs(list-value));
end
