function labels = subplot_labels(nn)
% Return a cell array of strings of the alphabet. 
% Useful for adding labels to subplots.
% If an integer is passed as an argument, then the label
% associated with that integer is returned.
%
% David Deepwell, 2018

labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)',...
          '(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)',...
          '(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};

if nargin == 1
    labels = labels{nn};
end
