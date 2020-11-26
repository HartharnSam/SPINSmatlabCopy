function T = SPINS_parameter_collator(list_filename)

%%
names_list = importdata(list_filename);
inds = strfind(names_list, ', ');
for i = 1:length(names_list)
    dir_list{i, 1} = names_list{i, 1}(1:inds{i, 1}(1)-1);
    slope_name{i,1} = names_list{i, 1}((inds{i, 1}(1)+2):inds{i, 1}(2)-1);
    amp_name{i,1} = names_list{i, 1}((inds{i, 1}(2)+2):inds{i, 1}(3)-1);
    category{i, 1} = names_list{i, 1}((inds{i, 1}(3)+2):end);
end

%%
loadonly = false; 
if exist('DATASummary.mat')
    load 'DATASummary.mat'
    if height(DATA) == length(dir_list)
        loadonly = true;
    end
end
if ~loadonly
cd('C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\Data\Model\NovakTank\020720_26');

params = spins_params;
chosen_params = [1:9 14:34 36:37 39:43];
load('wavestats.mat');

VarNames = fieldnames(params);
VarNames = [VarNames(chosen_params); fieldnames(WaveStats)];
VarTypes(1:length(VarNames), 1) = {'double'};

T = table('Size', [length(dir_list) length(VarNames)],'VariableTypes', VarTypes, 'VariableNames', VarNames, 'RowNames', dir_list);

for ii=1:length(dir_list)
    cd(['../', dir_list{ii}]);
    params = spins_params;
    load('wavestats.mat');
    
    for i = 1:length(VarNames)
        if any(i==1:length(chosen_params))
            paramvals(i) = extractfield(params, VarNames{i});
        else
            paramvals(i) = extractfield(WaveStats, VarNames{i});
        end
    end
    T(dir_list(ii), :) = num2cell(paramvals);
end

T = addvars(T, slope_name, 'Before', 'Lx');
T = addvars(T, amp_name, 'After', 'slope_name');
T = addvars(T, category, 'After', 'amp_name');

save('C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\Data\Model\NovakTank\DATASummary.mat', 'DATA');

else
    load 'DATASummary.mat'
    T = DATA;

end
