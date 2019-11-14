
clc

%% Hyperspectral unmixing

%%% Performs experiments from Section 4.1
disp('========================================================================')
disp('                        Hyperspectral unmixing                          ')
disp('========================================================================')

%%% set path to functions
clear all;
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:MKDIR:DirectoryExists')
rmpath('functions_geo_text')
addpath functions_hyperspectral
addpath(genpath('TOOLBOX_DWTRed_Frame'))

%%% set maximal duration and algorithmic method
time_max = 60;        % simulation duration in seconds
method   = 'PIPA-VM'; % 'PIPA' (without variable metric) or 'PIPA-VM'

%%% create result folder
folder   = 'results_hyperspectral'; % folder containing the results
filename = strcat(folder,'/','urban_',method); % file name
mkdir(folder);

%%% run PIPA or PIPA-VM
PIPA(method,time_max,filename);

%%% plot visual results, objective function and SNR with respect to time
my_plot(filename);

%% Joint geometry-texture decomposition and reconstruction of CT data

%%% Performs experiments from Section 4.2
disp('========================================================================')
disp('      Geometry-texture decomposition and reconstruction of CT data      ')
disp('========================================================================')

%%% set path to functions
clear all; 
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:MKDIR:DirectoryExists')
rmpath('functions_hyperspectral')
addpath functions_geo_text

%%% set maximal duration and sample name
sample_name = 'glass'; % sample name, 'glass' or 'agaricus'
time_max    = 5*60;    % simulation duration in seconds (about 30 min in the article) 

%%% create result folder
folder      ='results_geo_text'; % name of the folder containing the results
filename    = strcat(folder,'/',sample_name,'_PIPA-VM'); % file name
mkdir(folder);

%%% run PIPA-VM
PIPA(sample_name,time_max,filename);

%%% plot visual results, objective function and SNR with respect to time
my_plot(filename)
