function [] = par2var( params )
% Convert fields of the params structure into
% variables in the workspace that called this function

names = fieldnames(params);
for ii=1:length(names)
    val = params.(names{ii});
    assignin('caller', names{ii}, val)
end

