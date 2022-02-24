function [] = completion(ii, Nout, rate)
% completion prints the completed percentage of a script
%
% David Deepwell, 2019
% ii = Current frame
% Nout = Final frame
if nargin == 2
    rate = 0.1;
end

if ii == 1
    fprintf('Process 1 of %3d complete\n',Nout)
end

if rate < 1
    % print only a percentage of the times
    if mod(ii,Nout*rate) < 1
        do_print = true;
    else
        do_print = false;
    end
else
    % print every 'rate' steps
    if mod(ii,rate) == 0
        do_print = true;
    else
        do_print = false;
    end
end

if do_print
    fprintf('Process %3d of %3d complete: %3d%%\n',ii,Nout,round(ii/Nout*100))
end
