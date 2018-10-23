function check_make_dir(direc)
% CHECK_MAKE_DIR  Check if path exists in local directory
%                 and make the directory(ies) if it doesn't

if ~(exist([pwd,'/',direc],'dir')==7)
    mkdir(direc);
end
