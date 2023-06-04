clear all
clc
addpath([pwd,'/Functions'])
addpath([pwd,'/plate_models'])
addpath([pwd,'/input_model_normal'])
warning('off','all')

tStart = tic;
%% Load
model_file = 'input_model_30x30.mat';
load(model_file)
savename = 'vbplate_30x30_hbms';

%% Analysis setup
options.nModeI = 20;
options.nModeA = 30;
options.exppt = 12000;

%% Hierarchical reduction
uz_hbms = hss(uc_model, param, tree_model, options);

figure
freq = 0:2:1000;
semilogy(freq, abs(uz_hbms))

timing = toc(tStart);
save([savename,'.mat'],'freq','uz_hbms','timing')



