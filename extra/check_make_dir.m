function check_make_dir(the_path)
% CHECK_MAKE_DIR  Check if path exists in local directory
%                 and make the path if it doesn't exist

orig_dir = pwd();

% separate path into sub-directories
dirs = strsplit(the_path,'/');
% remove empty strings
dirs(strcmp('',dirs)) = [];

% loop over sub-directories
for nn = 1:length(dirs)
    direc = dirs{nn};
    if ~(exist(direc,'dir')==7)
        mkdir(direc);
    end
    cd(direc);
end

% return to original directory
cd(orig_dir);
