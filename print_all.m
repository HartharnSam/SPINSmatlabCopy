function [] = print_all(type, filenames)
% PRINT_ALL    Prints all the figures in the workspace
%
%  Usage:
%    print_figure()
%    print_figure('pdf')
%
%  Inputs:
%    'type'      - the file type
%    'filenames' - the name of the files
%
%  Outputs:
%    - none
%
% David Deepwepll, 2017

    % find all figures
    figs = findall(0,'type','figure');

    % check/make figures directory
    if ~exist('figures','dir');
        mkdir('figures')
    end

    % save to directory
    cd('figures')
    for ii = 1:length(figs)
        % choose figure and print
        fig_hand = figs(ii);
        figure(fig_hand)
        print_figure(fig_hand,[filenames,'_',num2str(fig_hand.Number)]);
        % print information
        disp(['Figure ',num2str(fig_hand.Number),' printed.'])
        completion(ii, length(figs))
    end
    cd('..')
end
