
fid = fopen([userpath, '\spinsstartup.m'], 'w');
fprintf(fid, '%s \n', ['current_path = ', pwd, ';']);
fprintf(fid, '%s \n', "%% Add SPINS_main/matlab to path");
fprintf(fid, '%s \n', "sm_matlab_loc = genpath(current_path);");
fprintf(fid, '%s \n', "addpath(sm_matlab_loc);");
fprintf(fid, '%s \n', "if ispc");
fprintf(fid, '%s \n', "    all_subs =  strsplit(sm_matlab_loc,';');");
fprintf(fid, '%s \n', "else");
fprintf(fid, '%s \n', "    all_subs = strsplit(sm_matlab_loc, ':');");
fprintf(fid, '%s \n', "end");
fprintf(fid, '%s \n', "for ii = 1:length(all_subs)");
fprintf(fid, '%s \n', '    if contains(all_subs{ii},"."+lettersPattern)');
fprintf(fid, '%s \n', '        rmpath(all_subs{ii});');
fprintf(fid, '%s \n', '    end');
fprintf(fid, '%s \n', 'end');
fclose(fid);

try
    open startup.m
catch
    fid = fopen([userpath, '\startup.m'], 'w');
    fclose(fid);
    open startup.m
end

disp('Add the following to the startup.m file, then close it:')
disp('addpath(genpath(usertoolbox))');
pause;

