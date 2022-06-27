clear all
close all
clc;
warning('off')

%%% ============ MAKE SURE THE PATHS ARE CORRECT =============== %%%
% Set path to the parant folder containing the analysis scripts downloaded from zenodo:
parentfolder = '/home/user07/Project/RippleDetection/';
% Set path for the required toolboxes:
path_to_toolboxes = '/home/user07/Project/RippleDetection/Toolbox';

% Add/remove paths: 
restoredefaultpath
%%% ============ MAKE SURE THE PATHS ARE CORRECT =============== %%%
addpath(fullfile(path_to_toolboxes,'fieldtrip-20220514'));
addpath(fullfile(parentfolder,'sub_routines_and_functions'));
addpath(genpath(fullfile(parentfolder,'general_code_and_data')));
set_figure_colors;