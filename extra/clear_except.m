function [] = clear_except(list)
% CLEAR_EXCEPT    Clear all variables except those in list
%
%  Usage:
%    clear_except({'a','b','c'})
%
%  Inputs:
%    'list' - a cell array containing strings of variables
%             that are to be kept
%
%  Outputs:
%    - none
%
% David Deepwell, 2017

caller_who = evalin('caller', 'who')';
clear_list = ['clear_list',setdiff(caller_who,list)];
assignin('caller', 'clear_list', clear_list)
evalin('caller', 'clear(clear_list{:})');
