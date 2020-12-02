%% spinsstartup script
% Sets up MATLAB environment to run SPINS stuff, execute spinsstartup to
% initiate
% Primarily, setting path for this session, secondly also producing a
% standard figure output "look"


%% add SPINSmatlab to path

sm_loc = genpath(['C:\Users\', getenv('username'), '\OneDrive - Newcastle University\Project\Shoal_Core\scripts\Model\SPINSmatlab']);
addpath(sm_loc);
if ispc
    all_subs =  strsplit(sm_loc,';');
else
    all_subs = strsplit(sm_loc, ':');
end
for ii = 1:length(all_subs)
    if ~isempty(strfind(all_subs{ii},'.'))
        rmpath(all_subs{ii});
    end
end

%% Add SPINS_main/matlab to path
sm_matlab_loc = genpath(['C:\Users\', getenv('username'), '\OneDrive - Newcastle University\Project\Shoal_Core\scripts\Model\SPINS_main\matlab']);
addpath(sm_matlab_loc);
if ispc
    all_subs =  strsplit(sm_matlab_loc,';');
else
    all_subs = strsplit(sm_matlab_loc, ':');
end
for ii = 1:length(all_subs)
    if ~isempty(strfind(all_subs{ii},'.'))
        rmpath(all_subs{ii});
    end
end


%% Set plotting variables as per "betterplots"
figure
betterplots(groot);

%% Clear workspace

cd(['C:\Users\', getenv('username'), '\OneDrive - Newcastle University\Project\Shoal_Core\Data\Model\']);
open SPINSContents
clearvars; close all; clc;

