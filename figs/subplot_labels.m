function labels = subplot_labels(nn, num_brackets)
% Return a cell array of strings of the alphabet. 
% Useful for adding labels to subplots.
% If an integer is passed as an argument, then the label
% associated with that integer is returned.
%
% David Deepwell, 2018

if ~exist('num_brackets', 'var')
    num_brackets = 2;
end

if ~isinteger(int8(num_brackets))
    error('Second argument must be the number of brackets')
end

if num_brackets == 1
    labels = {'a)','b)','c)','d)','e)','f)','g)','h)','i)',...
              'j)','k)','l)','m)','n)','o)','p)','q)','r)',...
              's)','t)','u)','v)','w)','x)','y)','z)'};
elseif num_brackets == 2
    labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)',...
              '(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)',...
              '(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};
end

if nargin > 0
    labels = labels{nn};
end
