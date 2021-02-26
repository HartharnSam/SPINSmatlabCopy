function dir_list = run_directory_names(LayerType, Selection)
%RUN_DIRECTORY_NAMES - Returns list of directories with SPINS runs in
%Runs are pre-entered in this function, to save long lists in many separate
%functions/scripts
%
%
% Inputs:
%    LayerType - REQUIRED: Type of stratification either "2 Layer", "3
%                       Layer" or "Continuous"
%    Selection - "full" returns the example subset for given stratification.
%                "Paper" returns the subset used in the JFM paper
%                "Examples" returns the subset used in fig 2 of JFM paper
%                   [OPTIONAL: Default = "full"]
%
% Outputs:
%    dir_list - List of runs as directory names for a given stratification
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 03-Dec-2020; Last revision: 03-Dec-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
if nargin<2
    Selection = "Full";
end

switch LayerType
    case "2 Layer"
        if strcmp(Selection, "Examples")
            dir_list = {'250720_31', '270520_04', '081020_44', '270720_33'};
        elseif strcmp(Selection, "Paper")
            dir_list = {'250720_31', '081020_44', '270720_33'};
        else
            dir_list = {'200520_02', '250520_03', '270520_04', '280520_05', '040620_06', ...
                '050620_07', '070620_08', '080620_09', '080620_10', '090620_11', '110620_12', ...
                '120620_13', '130620_14', '140620_15', '150620_16', '160620_17', '170620_18'...
                '230620_19', '240620_20', '250620_21', '260620_22', '270620_23', '280620_24',...
                '010720_25', '020720_26', '030720_27', '040720_28', '050720_29', '250720_31',...
                '260720_32', '270720_33', '300720_36', '081020_44'};
        end
    case "3 Layer"
        if strcmp(Selection, "Examples")
            dir_list = {'27_111120', '26_091120', '24_071020', '28_121120'};
        elseif strcmp(Selection, "Paper")
            dir_list = {'24_071020'};
        else
            dir_list = {'02_090720', '03_100720', '04_110720', '05_120720',...
                '07_310720', '24_071020', '25_221020', '26_091120', '27_111120',...
                '28_121120'};
        end
        
    case "Continuous"
        if strcmp(Selection, "Examples")
            dir_list = {'091020_45', '111020_47', '121020_48', '101120_49'};
        elseif strcmp(Selection, "Paper")
            dir_list = {'101020_46', '111020_47'};
        else % Full
            dir_list = {'091020_45', '101020_46', '111020_47', '121020_48',...
                '101120_49', '120220_50', '130220_51', '140220_52', '150220_53', ...
                '160220_54', '170220_55', '180220_56'};
        end
end
%---------------------------------------------------
%% END CODE %%
%---------------------------------------------------