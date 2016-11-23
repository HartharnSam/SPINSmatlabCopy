%cd /scratch/kglamb/ddeepwel/tankmode2/multidyes/2D/varypyc
cd /Users/ddeepwel/Orca_scratch/mode2shoal/Scotland/2D/w_hill_01/Lz035

%dirs = {'H006_dye','H007_dye','H008_dye','H010_dye','H012_dye'};
%dirs = {'H008_dye','H010_dye'};
dirs = {'H006','H008','H009','H010','H012','H014','H014_smalldt'};
%dirs = {'H008','H010','H012','H014'};

cd(dirs{1});
for mm = 1:length(dirs)
    clearvarlist = ['clearvarlist';setdiff(who,{'dirs';'mm'})];
    clear(clearvarlist{:}); 

    disp(['Directory: ',dirs{mm}])
    cd(['../',dirs{mm}])
    %diagnose
    %spins_tseries_multi
    spins_tseries_athill
end
