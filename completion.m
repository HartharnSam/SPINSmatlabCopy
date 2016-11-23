function [] = completion(ii, maxout, varargin)
% completion prints the completed percentage of a script
%
% David Deepwell, 2014

disp(['Process ',num2str(ii),' of ',...
      num2str(maxout),' complete, '...
      ,num2str(round(ii/maxout*100)),'%.'])
